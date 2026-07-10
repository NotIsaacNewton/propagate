//
// Created by Arian Dovald.
//

#include <string>
#include <numbers>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <chrono>
#include <print>
#include <algorithm>
#include "filetools.h"
#include "propagate.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"
#include "interpolate_1d.h"

// TODO: more safety checks and error paths (try <expected>)
//  allow for time-dependent potentials (write from potentials.cpp grid, add defineTDPotentialOperator)
//  implement CI models: propagation on multiple potential curves. make dimension-agnostic as much as possible.
//  generalize to higher dimensions

// creates array of squared momenta
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

// creates potential operator array from data file and outputs to op
void definePotentialOperator(const inputs& in, fftw_complex *op,
    const std::string& potfile, const bool imProp) {
    std::vector<double> potential(in.pot_grid);
    readArray1D(potfile, potential);
    // reshape potential with interpolation if needed
    if (in.space_grid != in.pot_grid) {
        std::vector<double> grid(in.pot_grid); // stores grid on which potential is defined
        const double dx = (in.final_pos-in.initial_pos)/(in.pot_grid-1); // width of potential grid
        // write potential grid
        for (int i = 0; i < in.pot_grid; i++) {
            grid[i] = in.initial_pos + i*dx;
        }
        spline_interp interpolator(grid, potential); // spline interpolation object
        std::vector<double> pot(in.space_grid); // temporary potential array
        // write temp array
        for (int i = 0; i < in.space_grid; i++) {
            pot[i] = interpolator.interp(in.initial_pos + i*in.dx);
        }
        // reshape potential
        potential = std::move(pot);
    }
    // write potential operator
    if (imProp) {
        // imaginary propagation
        for (int i = 0; i < in.space_grid; i++) {
            op[i][0] = exp(-potential[i] * in.dt / 2.0);
            op[i][1] = 0;
        }
    } else {
        // real propagation
        for (int i = 0; i < in.space_grid; i++) {
            const double phase = potential[i] * in.dt / 2.0;
            op[i][0] = cos(phase);
            op[i][1] = -sin(phase);
        }
    }
}

// calculates free-particle operator based on general values and outputs to op
void defineKineticOperator(const inputs& in, fftw_complex *op, const bool imProp) {
    const std::vector<double> mom = psquared(in.space_grid, in.dx);
    if (imProp) {
        // imaginary propagation
        for (int i = 0; i < in.space_grid; i++) {
            op[i][0] = exp(-in.dt * mom[i] / 2);
            op[i][1] = 0;
        }
    } else {
        // real propagation
        for (int i = 0; i < in.space_grid; i++) {
            const double phase = in.dt * mom[i] / 2;
            op[i][0] = cos(phase);
            op[i][1] = -sin(phase);
        }
    }
}

// applies potential energy exponential operator to psi
void applyPotentialOperator(const int gridpoints, fftw_complex *psi, const fftw_complex *V) {
    for (int i = 0; i < gridpoints; i++) {
        // sign flip for dft
        const int sign = i % 2 == 0 ? 1 : -1;
        psi[i][0] *= sign;
        psi[i][1] *= sign;
        // write
        const double re = psi[i][0];
        const double im = psi[i][1];
        // mult by V(x) exponential term
        psi[i][0] = re*V[i][0] - im*V[i][1];
        psi[i][1] = im*V[i][0] + re*V[i][1];
    }
}

// applies kinetic energy exponential operator to psi
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
    definePotentialOperator(in, V, potfile, imProp);
    defineKineticOperator(in, T, imProp);
    // return fftwResources
    return fftwResources{
        std::move(fft_ptr), std::move(ifft_ptr),
        std::unique_ptr<fftw_complex, void(*)(void*)>(V, fftw_free),
        std::unique_ptr<fftw_complex, void(*)(void*)>(T, fftw_free)
    };
}

// propagation tick
void propTick(const int gridpoints, fftw_complex *psi, const fftw_complex* V, const fftw_complex* T,
    fftw_plan fft, fftw_plan ifft, const double scale) {
    // apply e^(-i dt V(X))
    applyPotentialOperator(gridpoints, psi, V);
    // execute fft plan
    fftw_execute(fft);
    // apply e^(-i dt p^2 / 2)
    applyKineticOperator(gridpoints, psi, T);
    // execute inverse fft plan
    fftw_execute(ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    scale_fftw_complex(scale, psi, gridpoints);
    // apply e^(-i dt V(X))
    applyPotentialOperator(gridpoints, psi, V);
}

// propagates wavefunction based on general values
// Note: FFT methods implicitly impose periodic boundary conditions
void propagate(const inputs& in, fftw_complex *psi, const std::string& data, const bool imProp) {
    // prep fftw variables and plans
    auto [fft_ptr, ifft_ptr, Vp, Tp] = fftwPrep(in, psi, data, imProp);
    fftw_complex* V = Vp.get();
    fftw_complex* T = Tp.get();
    fftw_plan fft = fft_ptr.get();
    fftw_plan ifft = ifft_ptr.get();
    // scale for normalizing fft result
    const double scale = 1.0 / in.space_grid;
    // open output file
    const std::string output = data + "/psi_final.dat";
    std::ofstream wf(output, std::ios::app | std::ios::binary);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // prepare output buffer for entire set of points
    std::vector<double> buffer;
    buffer.reserve((in.time_grid / in.nt_prints + 1) * (in.space_grid / in.nx_prints) * 2);
    // norm psi (often slightly off norm) and check norm
    std::vector<double> psi_squared(in.space_grid); // stores |psi|^2
    fftw_complex_square(psi, psi_squared); // calculates |psi|^2
    double norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared); // calculate and store norm
    scale_fftw_complex(1/sqrt(norm), psi, in.space_grid); // normalize psi
    fftw_complex_square(psi, psi_squared); // recalculate |psi|^2
    norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared); // recalculate norm
    std::cout << RED << "Check norm:\n" << RESET;
    std::print("The initial norm is {}\n", norm);
    // spacer
    spacer(RESET);
    // console output
    std::cout << "Propagation progress:\n";
    // propagation loop
    for (int t = 0; t <= in.time_grid; t++) {
        // print completion % to console
        !(t % (in.time_grid / 10)) ? std::cout << GREEN << "\r" << 100*t/in.time_grid << "%" : std::cout << RESET;
        // write lines in output file
        writeOutput(psi, t, in.space_grid, in.nx_prints, in.nt_prints, buffer);
        // propagate for one tick
        propTick(in.space_grid, psi, V, T, fft, ifft, scale);
        // naive renormalization if doing imaginary time propagation
        if (imProp && t % 10 == 0) {
            fftw_complex_square(psi, psi_squared);
            norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared);
            scale_fftw_complex(1/sqrt(norm), psi, in.space_grid);
        }
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

    // warning if odd number of gridpoints
    if (in.space_grid % 2 != 0) {
        spacerFancy(YELLOW);
        std::cerr << YELLOW << "WARNING: Odd number of gridpoints may produce asymmetric momentum grid.\n" << RESET;
        spacerFancy(YELLOW);
    }
    // warning if large dt
    if (in.time_grid < 1000) {
        spacerFancy(YELLOW);
        std::cerr << YELLOW << "WARNING: Large dt may fail to resolve rapid dynamics. Use dt < pi/E_max.\n" << RESET;
        spacerFancy(YELLOW);
    }

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
