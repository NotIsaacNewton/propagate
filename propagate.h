//
// Created by Arian Dovald on 6/23/25.
//

#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <functional>
#include <numbers>
#include "filetools.h"
#include "fftw3.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"

// calculates momentum-squared based on general values
static double psquared(int i, int gridpoints, double space_width) {
    double p = 2 * std::numbers::pi * i / (gridpoints * space_width) - (std::numbers::pi / space_width);
    return p*p;
}

// real time propagation
// calculates free-particle operator based on general values and outputs to op
static void defineKineticOperator(int gridpoints, fftw_complex *op, double space_width, double time_width) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = cos(time_width * psquared(i, gridpoints, space_width) / 2);
        op[i][1] = -sin(time_width * psquared(i, gridpoints, space_width) / 2);
    }
}

// potential energy exponential operator
static void potentialOperator(int gridpoints, fftw_complex *psi, const fftw_complex *V) {
    for (int i = 0; i < gridpoints; i++) {
        int sign = (i % 2 == 0) ? 1 : -1;
        // sign change for dft
        psi[i][0] *= sign;
        psi[i][1] *= sign;
        double re = psi[i][0];
        double im = psi[i][1];
        // mult by V(x) exponential term
        psi[i][0] = re*V[i][0] - im*V[i][1];
        psi[i][1] = im*V[i][0] + re*V[i][1];
    }
}

// kinetic energy exponential operator
static void kineticOperator(int gridpoints, fftw_complex *psi, const fftw_complex *T) {
    for (int i = 0; i < gridpoints; i++) {
        double re = psi[i][0];
        double im = psi[i][1];
        psi[i][0] = re * T[i][0] - im * T[i][1];
        psi[i][1] = im * T[i][0] + re * T[i][1];
    }
}

// output to file
static void writeOutput(const fftw_complex *psi, int t, double start, int gridpoints, double space_width, double time_width,
                 int output_ndx, int output_ndt, std::ostringstream &buffer) {
    // print to output file every output_ndt steps
    if (t % output_ndt == 0) {
        // print to ouput file every output_ndx steps
        for (int i=0; i<gridpoints; i += output_ndx) {
            // prepare output
            double re = psi[i][0];
            double im = psi[i][1];
            double square = re*re + im*im;
            // send stuff to buffer
            buffer << t * time_width << " " << i*space_width + start << " " << re << " " << im << " " << square << "\n";
        }
    }
}

// propagates wavefunction based on general values
// TODO:
//  allow for time-dependent potentials (put in definePotentialOperator grid)
//  find some way to dynamically define re and im, to reduce redundancy
[[maybe_unused]] inline void propagate(const inputs& in,
    fftw_complex *psi, const std::string& data) {
    // create fft and inverse fft plans with RAII
    std::string wisdomfile = data + "/fftw_wisdom.dat";
    fftw_import_wisdom_from_filename(wisdomfile.c_str());
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> fft_ptr(
            fftw_plan_dft_1d(in.space_grid, psi, psi, FFTW_FORWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> ifft_ptr(
            fftw_plan_dft_1d(in.space_grid, psi, psi, FFTW_BACKWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    fftw_export_wisdom_to_filename(wisdomfile.c_str());
    // define potential term
    std::string potfile = data + "/potential.dat";
    auto V = fftw_alloc_complex(in.space_grid);
    fftw_complex_array_from_file(potfile, V);
    // define kinetic term;
    auto T = fftw_alloc_complex(in.space_grid);
    defineKineticOperator(in.space_grid, T, in.dx, in.dt);
    // RAII wrapper for T and V
    std::unique_ptr<fftw_complex, void(*)(void*)> Vp{V, fftw_free}, Tp{T, fftw_free};
    // scale for normalizing fft result
    double scale = 1.0 / in.space_grid;
    // open output file
    std::string output = data + "/psi_final.dat";
    std::ofstream wf;
    wf.open(output, std::ios::app);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // prepare input buffer for entire set of points
    std::ostringstream buffer;
    // check norm
    double psi_squared[in.space_grid];
    fftw_complex_square(in.space_grid, psi, psi_squared);
    std::cout << RED << "Check norm:\n" << RESET;
    double norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared);
    std::print("The initial norm is {}\n", norm);
    // spacer
    spacer(RESET);
    // console output
    std::cout << "Propagation progress:\n";
    // propagation loop
    for (int t = 0; t <= in.time_grid; t++) {
        // print completion % to console
        !(t % static_cast<int>(in.time_grid * 0.1)) ?
        std::cout << GREEN << "\r" << 100*t/in.time_grid << "%" : std::cout << RESET;
        // write lines in output file
        writeOutput(psi, t, in.initial_pos, in.space_grid, in.dx, in.dt,
            in.nx_prints, in.nt_prints, buffer);
        // apply e^(-i dt V(X))
        potentialOperator(in.space_grid, psi, V);
        // execute fft plan
        fftw_execute(fft_ptr.get());
        // apply e^(-i dt p^2 / 2)
        kineticOperator(in.space_grid, psi, T);
        // execute inverse fft plan
        fftw_execute(ifft_ptr.get());
        // normalize fftw result (fftw uses non-normalized fft algorithm)
        scale_fftw_complex(scale, psi, in.space_grid);
        // apply e^(-i dt V(X))
        potentialOperator(in.space_grid, psi, V);
        // execute fft plan
    }
    // save buffer to output file and close the file
    wf << buffer.str();
    wf.close();
    // console output
    std::cout << "\n";
    // spacer
    spacer(RESET);
    // check norm
    fftw_complex_square(in.space_grid, psi, psi_squared);
    std::cout << RED << "Check norm:\n" << RESET;
    norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared);
    std::print("The final norm is {}\n", norm);
}

// imaginary time propagation
// calculates free-particle operator based on general values and outputs to op
static void defineImKineticOperator(int gridpoints, fftw_complex *op, double space_width, double time_width) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = exp(-time_width * psquared(i, gridpoints, space_width) / 2);
        op[i][1] = 0;
    }
}

// creates potential operator array from potential function and outputs to op
static void defineImPotentialOperator(int gridpoints, fftw_complex *op, const std::function<double(double)>& potential,
                               double space_width, double time_width, double start) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = exp(-potential(i * space_width + start) * time_width / 2.0);
        op[i][1] = 0;
    }
}

// propagates wavefunction based on general values
[[maybe_unused]] inline void ipropagate(const inputs& in,
    fftw_complex *psi, const std::function<double(double)>& potential, std::string& data) {
    // create fft and inverse fft plans with RAII
    std::string wisdomfile = data + "/fftw_wisdom.dat";
    fftw_import_wisdom_from_filename(wisdomfile.c_str());
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> fft_ptr(
            fftw_plan_dft_1d(in.time_grid, psi, psi, FFTW_FORWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> ifft_ptr(
            fftw_plan_dft_1d(in.time_grid, psi, psi, FFTW_BACKWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    fftw_export_wisdom_to_filename(wisdomfile.c_str());
    // define potential term
    auto V = fftw_alloc_complex(in.time_grid);
    defineImPotentialOperator(in.time_grid, V, potential, in.dx, in.dt, in.initial_pos);
    // define e^(-i*dt*T) operator;
    auto T = fftw_alloc_complex(in.time_grid);
    defineImKineticOperator(in.time_grid, T, in.dx, in.dt);
    // RAII wrapper for T and V
    std::unique_ptr<fftw_complex, void(*)(void*)> Vp{V, fftw_free}, Tp{T, fftw_free};
    // scale for normalizing fft results
    double scale = 1.0 / in.time_grid;
    // open output file
    std::string output = data + "/psi_final.dat";
    std::ofstream wf;
    wf.open(output, std::ios::app);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // prepare input buffer for entire set of points
    std::ostringstream buffer;
    // console output
    std::cout << "Propagation progress:\n";
    // propagation loop
    for (int t = 0; t <= in.time_grid; t++) {
        // print completion % to console
        !(t % static_cast<int>(in.time_grid * 0.1)) ?
        std::cout << GREEN << "\r" << 100*t/in.time_grid << "%" : std::cout << RESET;
        // write lines in output file
        writeOutput(psi, t, in.initial_pos, in.time_grid, in.dx, in.dt,
            in.nx_prints, in.nt_prints, buffer);
        // apply e^(-i dt V(X))
        potentialOperator(in.time_grid, psi, V);
        // execute fft plan
        fftw_execute(fft_ptr.get());
        // apply e^(-i dt p^2 / 2)
        kineticOperator(in.time_grid, psi, T);
        // execute inverse fft plan
        fftw_execute(ifft_ptr.get());
        // normalize fftw result (fftw uses non-normalized fft algorithm)
        scale_fftw_complex(scale, psi, in.time_grid);
        // apply e^(-i dt V(X))
        potentialOperator(in.time_grid, psi, V);
    }
    // save buffer to output file and close the file
    wf << buffer.str();
    wf.close();
    // console output
    std::cout << "\n";
}

#endif //PROPAGATE_H
