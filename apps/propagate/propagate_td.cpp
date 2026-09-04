//
// Created by Arian Dovald on 8/19/26.
//

#include "propagate_td.h"
#include "console_tools.h"
#include "fftw_complex_tools.h"
#include <fstream>
#include <memory>
#include <print>
#include <chrono>
#include <iostream>

// gets potential and returns arrays
std::vector<std::vector<double>> getPotential(const inputs& in, const std::string& data) {
    // prepare potential arrays
    std::vector<std::vector<double>> potential; // full potential array
    readArray2D(data + "/potential.dat", potential,
        in.space_grid, in.time_grid); // reads in potential
    return potential;
}

// creates potential operator array from data and outputs to op
void definePotentialOperatorTD(const inputs& in, const int& tick, fftw_complex *op,
    const std::vector<std::vector<double>>& potential) {
    // write potential operator
    for (int i = 0; i < in.space_grid; i++) {
        const double phase = potential[tick][i] * in.dt / 2.0;
        op[i][0] = cos(phase);
        op[i][1] = -sin(phase);
    }
}

// prepares fftw variables with RAII
fftwResources fftwPrepTD(const inputs& in, fftw_complex *psi, const std::string& data) {
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
    // save fftw wisdom to file
    fftw_export_wisdom_to_filename(wisdomfile.c_str());
    // locate potential curve file
    const std::string potfile = data + "/potential.dat";
    // define only kinetic term since potential is defined at each time step
    const auto V = fftw_alloc_complex(in.space_grid);
    const auto T = fftw_alloc_complex(in.space_grid);
    defineKineticOperator(in, T, false);
    // return fftwResources
    return fftwResources{
        .fft_ptr = std::move(fft_ptr), .ifft_ptr = std::move(ifft_ptr),
        .Vp = std::unique_ptr<fftw_complex, void(*)(void*)>(V, fftw_free),
        .Tp = std::unique_ptr<fftw_complex, void(*)(void*)>(T, fftw_free)
    };
}

// propagation tick
void propTickTD(const int& tick, const inputs& in, fftw_complex *psi, const std::vector<std::vector<double>>& potential,
    fftw_complex* V, const fftw_complex* T, const fftw_plan* fft, const fftw_plan* ifft, const double scale) {
    // define potential operator
    definePotentialOperatorTD(in, tick, V, potential);
    // apply e^(-i dt V(x,t))
    applyPotentialOperator(in.space_grid, psi, V);
    // execute fft plan
    fftw_execute(*fft);
    // apply e^(-i dt p^2 / 2)
    applyKineticOperator(in.space_grid, psi, T);
    // execute inverse fft plan
    fftw_execute(*ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    scale_fftw_complex(scale, psi, in.space_grid);
    // apply e^(-i dt V(x))
    applyPotentialOperator(in.space_grid, psi, V);
}

// propagates wavefunction based on general values
void propagateTD(const inputs& in, const std::string& data) {
    // get initial wavepacket
    auto psip = getWavepacket(in, data);
    auto psi = psip.get();
    // prep fftw variables and plans
    auto [fft_ptr, ifft_ptr, Vp, Tp] = fftwPrepTD(in, psi, data);
    fftw_complex* V = Vp.get();
    fftw_complex* T = Tp.get();
    fftw_plan fft = fft_ptr.get();
    fftw_plan ifft = ifft_ptr.get();
    // get potential
    auto potential = getPotential(in, data);
    // scale for normalizing fft result
    const double scale = 1.0 / in.space_grid;
    // open output file and buffer
    auto [wf, buffer] = openWFOutputFile(in, data);
    // norm psi (often slightly off norm) and check norm
    normalize(in.space_grid, in.dx, psi); // normalize psi
    double mag = norm(in.space_grid, in.dx, psi); // recalculate norm
    std::cout << RED << "Check norm:\n" << RESET;
    std::print("The initial norm is {}\n", mag);
    // spacer
    spacer(RESET);
    // console output
    std::cout << "Propagating...\n";
    // propagation loop
    for (int t = 0; t < in.time_grid; t++) {
        // print completion % to console
        !((t+1) % (in.time_grid / 10)) ? progressBar(GREEN, 100*(t+1)/in.time_grid) : reset();
        // write lines in output file
        writeOutput(psi, t, in.space_grid, in.nx_prints, in.nt_prints, buffer);
        // propagate for one tick
        propTickTD(t, in, psi, potential, V, T, &fft, &ifft, scale);
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
    mag = norm(in.space_grid, in.dx, psi);
    std::cout << RED << "Check norm:\n" << RESET;
    std::print("The final norm is {}\n", mag);
}
