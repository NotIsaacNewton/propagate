//
// Created by Arian Dovald.
//

#include <string>
#include <sstream>
#include <numbers>
#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <memory>
#include <chrono>
#include <print>
#include "filetools.h"
#include "propagate.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"

// TODO: more safety checks and error paths (try <expected>)
//  interpolation routines that fit any potential to the grid parameters
//  allow for time-dependent potentials (write from potentials.cpp grid, add defineTDPotentialOperator)
//  block buffers
//  generalize to higher dimensions

std::vector<double> psquared(const int gridpoints, const double space_width) {
    const double scale = 2 * std::numbers::pi / (gridpoints * space_width);
    const double shift = - std::numbers::pi / space_width;
    std::vector<double> mom(gridpoints);
    for (int i = 0; i < gridpoints; i++) {
        mom[i] = scale * i + shift;
        mom[i] *= mom[i];
    }
    return mom;
}

// real time propagation
// creates potential operator array from data file and outputs to op
void definePotentialOperator(const int gridpoints, fftw_complex *op,
    const std::string& potfile, const double time_width, const bool imProp) {
    std::vector<double> potential(gridpoints);
    readArray1D(potfile, potential);
    if (imProp) {
        for (int i = 0; i < gridpoints; i++) {
            op[i][0] = exp(-potential[i] * time_width / 2.0);
            op[i][1] = 0;
        }
    } else {
        for (int i = 0; i < gridpoints; i++) {
            const double phase = potential[i] * time_width / 2.0;
            op[i][0] = cos(phase);
            op[i][1] = -sin(phase);
        }
    }
}

// calculates free-particle operator based on general values and outputs to op
void defineKineticOperator(const int gridpoints, fftw_complex *op, const double space_width,
    const double time_width, const bool imProp) {
    const std::vector<double> mom = psquared(gridpoints, space_width);
    if (imProp) {
        for (int i = 0; i < gridpoints; i++) {
            op[i][0] = exp(-time_width * mom[i] / 2);
            op[i][1] = 0;
        }
    } else {
        for (int i = 0; i < gridpoints; i++) {
            const double phase = time_width * mom[i] / 2;
            op[i][0] = cos(phase);
            op[i][1] = -sin(phase);
        }
    }
}

// potential energy exponential operator
void applyPotentialOperator(const int gridpoints, fftw_complex *psi, const fftw_complex *V) {
    for (int i = 0; i < gridpoints; i++) {
        const int sign = i % 2 == 0 ? 1 : -1;
        // sign change for dft
        psi[i][0] *= sign;
        psi[i][1] *= sign;
        const double re = psi[i][0];
        const double im = psi[i][1];
        // mult by V(x) exponential term
        psi[i][0] = re*V[i][0] - im*V[i][1];
        psi[i][1] = im*V[i][0] + re*V[i][1];
    }
}

// kinetic energy exponential operator
void applyKineticOperator(const int gridpoints, fftw_complex *psi, const fftw_complex *T) {
    for (int i = 0; i < gridpoints; i++) {
        const double re = psi[i][0];
        const double im = psi[i][1];
        psi[i][0] = re * T[i][0] - im * T[i][1];
        psi[i][1] = im * T[i][0] + re * T[i][1];
    }
}

// output to file
void writeOutput(const fftw_complex *psi, const int t, const int gridpoints, const int output_ndx, const int output_ndt,
    std::vector<double> &buffer) {
    // print to output file every output_ndt steps
    if (t % output_ndt == 0) {
        // print to output file every output_ndx steps
        for (int i=0; i<gridpoints; i += output_ndx) {
            // send stuff to buffer
            buffer.push_back(psi[i][0]);
            buffer.push_back(psi[i][1]);
        }
    }
}

// struct bundling all resources needed for using fftw in propagation loops
struct fftwResources {
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> fft_ptr;
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> ifft_ptr;
    std::unique_ptr<fftw_complex, void(*)(void*)> Vp;
    std::unique_ptr<fftw_complex, void(*)(void*)> Tp;
};

// prepares fftw variables with RAII
fftwResources fftwPrep(const inputs& in, fftw_complex *psi, const std::string& data, const bool imProp) {
    // locate fftw wisdom file
    const std::string wisdomfile = data + "/fftw_wisdom.dat";
    fftw_import_wisdom_from_filename(wisdomfile.c_str());
    // create fft and inverse fft plans with RAII
    auto fft_ptr = std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>(
            fftw_plan_dft_1d(in.space_grid, psi, psi, FFTW_FORWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    auto ifft_ptr = std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>(
            fftw_plan_dft_1d(in.space_grid, psi, psi, FFTW_BACKWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    // save fftw wisom to file
    fftw_export_wisdom_to_filename(wisdomfile.c_str());
    // locate potential curve file
    const std::string potfile = data + "/potential.dat";
    // define potential and kinetic terms
    const auto V = fftw_alloc_complex(in.space_grid);
    const auto T = fftw_alloc_complex(in.space_grid);
    definePotentialOperator(in.space_grid, V, potfile, in.dt, imProp);
    defineKineticOperator(in.space_grid, T, in.dx, in.dt, imProp);
    // return fftwResources
    return fftwResources{
        std::move(fft_ptr), std::move(ifft_ptr),
        std::unique_ptr<fftw_complex, void(*)(void*)>(V, fftw_free),
        std::unique_ptr<fftw_complex, void(*)(void*)>(T, fftw_free)
    };
}

// propagation tick
void propTick(const inputs& in, fftw_complex *psi, const fftw_complex* V, const fftw_complex* T,
    fftw_plan fft, fftw_plan ifft, const double scale) {
    // apply e^(-i dt V(X))
    applyPotentialOperator(in.space_grid, psi, V);
    // execute fft plan
    fftw_execute(fft);
    // apply e^(-i dt p^2 / 2)
    applyKineticOperator(in.space_grid, psi, T);
    // execute inverse fft plan
    fftw_execute(ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    scale_fftw_complex(scale, psi, in.space_grid);
    // apply e^(-i dt V(X))
    applyPotentialOperator(in.space_grid, psi, V);
}

// propagates wavefunction based on general values
void propagate(const inputs& in, fftw_complex *psi, const std::string& data, bool imProp) {
    // prep fftw variables and plans
    auto [fft_ptr, ifft_ptr, Vp, Tp] = fftwPrep(in, psi, data, imProp);
    fftw_complex* V = Vp.get();
    fftw_complex* T = Tp.get();
    fftw_plan fft = fft_ptr.get();
    fftw_plan ifft = ifft_ptr.get();
    // scale for normalizing fft result
    const double scale = 1.0 / in.space_grid;
    // open output file
    std::string output = data + "/psi_final.dat";
    std::ofstream wf;
    wf.open(output, std::ios::app);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // prepare input buffer for entire set of points
    std::vector<double> buffer;
    buffer.reserve((in.time_grid / in.nt_prints + 1) * (in.space_grid / in.nx_prints) * 2);
    // check norm
    std::vector<double> psi_squared(in.space_grid);
    fftw_complex_square(psi, psi_squared);
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
        writeOutput(psi, t, in.space_grid, in.nx_prints, in.nt_prints, buffer);
        // propagate for one tick
        propTick(in, psi, V, T, fft, ifft, scale);
    }
    // save buffer to output file and close the file
    wf.write(reinterpret_cast<const char*>(buffer.data()),
         static_cast<std::streamsize>(buffer.size() * sizeof(double)));
    wf.close();
    // console output
    std::cout << "\n";
    // spacer
    spacer(RESET);
    // check norm
    fftw_complex_square(psi, psi_squared);
    std::cout << RED << "Check norm:\n" << RESET;
    norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared);
    std::print("The final norm is {}\n", norm);
}

// inputs: location/of/input_file location/of/data_directory
int main(const int argc, const char* argv[]) {
    if (argc != 4) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}[location/of/input_file] [location/of/data_directory] [im prop: true/false] \n", GREEN);
        spacerFancy(RED);
        return 1;
    }

    // record the start time
    auto start = std::chrono::steady_clock::now();

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string psifile = data + "/psi_initial.dat";
    const std::string psiout = data + "/psi_final.dat";

    // imaginary propagation?
    const bool imProp = std::string(argv[3]) == "true";

    // spacer
    spacerChunky(BLUE);
    std::print("\n");

    // read input file
    const inputs in = readInputs(inputfile);

    // clear output file
    std::ofstream clearout;
    clearout.open(psiout);
    clearout.close();

    // allocate wavefunction with RAII
    auto psi = fftw_alloc_complex(in.space_grid);
    std::unique_ptr<fftw_complex, void(*)(void*)> psip{psi, fftw_free};
    // get wavefunction from psifile and save it to psi
    fftw_complex_array_from_file(psifile, psi, in.space_grid);

    // spacer
    spacerThick(RESET);

    // propagate wf
    propagate(in, psi, data, imProp);

    // spacer
    spacerThick(RESET);

    // record end time and duration
    auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    std::print("Done :)\n\n");
    std::print("{}Execution time: {} seconds\n", BLUE, sec.count());

    // spacer
    spacerChunky(BLUE);

    return 0;
}
