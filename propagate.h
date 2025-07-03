//
// Created by Arian Dovald on 6/23/25.
//

#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <functional>
#include "armpl.h"
#include "fftw3.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"

// scales entire array by a scalar
void scale_fftw_complex(double scalar, fftw_complex *complex_vec, int size) {
    for (int i = 0; i < size; i++) {
        complex_vec[i][0] *= scalar;
        complex_vec[i][1] *= scalar;
    }
}

// calculates momentum-squared based on general values
double psquared(int i, int gridpoints, double space_width) {
    double p = 2 * M_PI * i / (gridpoints * space_width) - (M_PI / space_width);
    return p*p;
}

// real time propagation
// calculates free-particle operator based on general values and outputs to op
void defineKineticOperator(int gridpoints, fftw_complex *op, double space_width, double time_width) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = cos(time_width * psquared(i, gridpoints, space_width) / 2);
        op[i][1] = -sin(time_width * psquared(i, gridpoints, space_width) / 2);
    }
}

// creates potential operator array from potential function and outputs to op
void definePotentialOperator(int gridpoints, fftw_complex *op, const std::function<double(double)>& potential,
                             double space_width, double time_width, double start) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = cos(potential(i * space_width + start) * time_width / 2.0);
        op[i][1] = -sin(potential(i * space_width + start) * time_width / 2.0);
    }
}

// potential energy exponential operator
void potentialOperator(int gridpoints, fftw_complex *psi, fftw_complex *V) {
    int sign;
    for (int i = 0; i < gridpoints; i++) {
        sign = (i % 2 == 0) ? 1 : -1;
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
void kineticOperator(int gridpoints, fftw_complex *psi, fftw_complex *T) {
    for (int i = 0; i < gridpoints; i++) {
        double re = psi[i][0];
        double im = psi[i][1];
        psi[i][0] = re * T[i][0] - im * T[i][1];
        psi[i][1] = im * T[i][0] + re * T[i][1];
    }
}

// output to file
void writeOutput(std::ofstream& wf, fftw_complex* psi, int t, double start, int gridpoints,
            double space_width, double time_width, int output_ndx, int output_ndt) {
    // print to output file every t_out steps
    if (t % output_ndt == 0) {
        // prepare input buffer for entire set of points
        std::ostringstream  buffer;
        for (int i=0; i<gridpoints; i += output_ndx) {
            // prepare output
            double re = psi[i][0];
            double im = psi[i][1];
            double square = re*re + im*im;
            // send stuff to buffer
            buffer << t * time_width << " " << i*space_width + start << " " << re << " " << im << " " << square << "\n";
        }
        // save buffer to file
        wf << buffer.str();
    }
}

// propagates wavefunction based on general values
// TODO:
//  collect inputs into struct or something (this is becoming clunky and annoying to work with)
//  allow for time-dependent potentials (put in definePotentialOperator grid)
//  use multi-threaded fftw methods :)
//  find some way to dynamically define re and im, to reduce redundancy
//  wrap plan creation and destruction in some kind of clean helper function —
//  FIXME: need RAII
[[maybe_unused]] void propagate(int gridpoints, double space_width, double start, fftw_complex *psi,
                                const std::function<double(double)>& potential, double time_width, int steps,
                                std::string& output, int output_ndt, int output_ndx) {
    // create fft and inverse fft plans
    fftw_plan fft = fftw_plan_dft_1d(gridpoints, psi, psi, -1, FFTW_ESTIMATE);
    fftw_plan ifft = fftw_plan_dft_1d(gridpoints, psi, psi, 1, FFTW_ESTIMATE);
    // define potential term (once because it is time-independent)
    fftw_complex V[gridpoints];
    definePotentialOperator(gridpoints, V, potential, space_width, time_width, start);
    // define kinetic term (once because it is time-independent);
    fftw_complex T[gridpoints];
    defineKineticOperator(gridpoints, T, space_width, time_width);
    // scale for normalizing fft result
    double scale = 1.0 / gridpoints;
    // open output file
    std::ofstream wf;
    wf.open(output, std::ios::app);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // check norm
    double psi_squared[gridpoints];
    fftw_complex_square(gridpoints, psi, psi_squared);
    std::cout << RED << "Check norm:\n" << RESET;
    double norm = fftw_complex_integrate(gridpoints, space_width, psi_squared);
    std::print("The initial norm is {}\n", norm);
    // spacer
    spacer(RESET);
    // console output
    std::cout << "Propagation progress:\n";
    // propagation loop
    for (int t = 0; t <= steps; t++) {
        // print completion % to console
        std::cout << GREEN << "\r" << 100*t/steps << "%";
        // write lines in output file
        writeOutput(wf, psi, t, start, gridpoints, space_width, time_width, output_ndx, output_ndt);
        // apply e^(-i dt V(X))
        potentialOperator(gridpoints, psi, V);
        // execute fft plan
        fftw_execute(fft);
        // apply e^(-i dt p^2 / 2)
        kineticOperator(gridpoints, psi, T);
        // execute inverse fft plan
        fftw_execute(ifft);
        // normalize fftw result (fftw uses non-normalized fft algorithm)
        scale_fftw_complex(scale, psi, gridpoints);
        // apply e^(-i dt V(X))
        potentialOperator(gridpoints, psi, V);
        // execute fft plan
    }
    // close output file
    wf.close();
    // console output
    std::cout << "\n";
    // clean up
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    // spacer
    spacer(RESET);
    // check norm
    fftw_complex_square(gridpoints, psi, psi_squared);
    std::cout << RED << "Check norm:\n" << RESET;
    norm = fftw_complex_integrate(gridpoints, space_width, psi_squared);
    std::print("The final norm is {}\n", norm);
}

// imaginary time propagation
// calculates free-particle operator based on general values and outputs to op
void defineImKineticOperator(int gridpoints, fftw_complex *op, double space_width, double time_width) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = exp(-time_width * psquared(i, gridpoints, space_width) / 2);
        op[i][1] = 0;
    }
}

// creates potential operator array from potential function and outputs to op
void defineImPotentialOperator(int gridpoints, fftw_complex *op, const std::function<double(double)>& potential,
                               double space_width, double time_width, double start) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = exp(-potential(i * space_width + start) * time_width / 2.0);
        op[i][1] = 0;
    }
}

// propagates wavefunction based on general values
[[maybe_unused]] void ipropagate(int gridpoints, double space_width, double start, fftw_complex *psi,
                                 const std::function<double(double)>& potential, double time_width, int steps,
                                 std::string& output, int output_ndt, int output_ndx) {
    // create fft and inverse fft plans
    fftw_plan fft = fftw_plan_dft_1d(gridpoints, psi, psi, -1, FFTW_ESTIMATE);
    fftw_plan ifft = fftw_plan_dft_1d(gridpoints, psi, psi, 1, FFTW_ESTIMATE);
    // define potential term (once because it is time-independent)
    fftw_complex V[gridpoints];
    defineImPotentialOperator(gridpoints, V, potential, space_width, time_width, start);
    // define e^(-i*dt*T) operator (once because it is time-independent);
    fftw_complex T[gridpoints];
    defineImKineticOperator(gridpoints, T, space_width, time_width);
    // scale for normalizing fft results
    double scale = 1.0 / gridpoints;
    // open output file
    std::ofstream wf;
    wf.open(output, std::ios::app);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // console output
    std::cout << "Propagation progress:\n";
    // propagation loop
    for (int t = 0; t <= steps; t++) {
        // print completion % to console
        std::cout << GREEN << "\r" << 100*t/steps << "%";
        // write lines in output file
        writeOutput(wf, psi, t, start, gridpoints, space_width, time_width, output_ndx, output_ndt);
        // apply e^(-i dt V(X))
        potentialOperator(gridpoints, psi, V);
        // execute fft plan
        fftw_execute(fft);
        // apply e^(-i dt p^2 / 2)
        kineticOperator(gridpoints, psi, T);
        // execute inverse fft plan
        fftw_execute(ifft);
        // normalize fftw result (fftw uses non-normalized fft algorithm)
        scale_fftw_complex(scale, psi, gridpoints);
        // apply e^(-i dt V(X))
        potentialOperator(gridpoints, psi, V);
    }
    // close output file
    wf.close();
    // console output
    std::cout << "\n";
    // clean up
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
}

#endif //PROPAGATE_H
