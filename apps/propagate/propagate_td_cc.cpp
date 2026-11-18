//
// Created by Arian Dovald on 9/2/26.
//

#include "propagate_td_cc.h"
#include "fftw_complex_tools.h"
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <print>
#include <algorithm>
#include <variant>
#include "console_tools.h"

// gets potential and returns arrays
std::vector<std::vector<hermitian_matrix> > getPotentialCC(const inputs &in, const std::string &data) {
    // prepare potential arrays
    std::vector potential(in.time_grid,
    std::vector(in.space_grid, hermitian_matrix(in.channels))); // full potential array
    std::vector<std::vector<double>> temp; // array containing raw values
    std::print("Reading potential from {}/potential.dat...\n", data);
    readArray2D(data + "/" + in.pot_file, temp,
                in.space_grid*(in.channels*in.channels+in.channels),
                in.time_grid); // reads in potential file
    #pragma omp parallel for default(none) shared(in, potential, temp)
    for (int i = 0; i < in.time_grid; i++) {
        for (int c1 = 0; c1 < in.channels; c1++) {
            for (int c2 = c1; c2 < in.channels; c2++) {
                for (int j = 0; j < in.space_grid; j++) {
                    potential[i][j](c1,c2) = std::complex(
                        temp[i][2*(j + c1*(in.channels-(c1+1)/2)*in.space_grid + c2*in.space_grid)],
                        temp[i][2*(j + c1*(in.channels-(c1+1)/2)*in.space_grid + c2*in.space_grid)+1]
                        );
                }
            }
        }
    }
    std::print("\nPotential read!\n");
    return potential;
}

// creates array of diagonal potential operator arrays from potential at tick and outputs to op
void definePotentialOperatorDiag(const inputs& in, const int& tick, const std::vector<fftw_complex*>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential) {
    #pragma omp parallel for default(none) shared(in, potential, tick, op)
    for (int i = 0; i < in.space_grid; i++) {
        for (int c = 0; c < in.channels; c++) {
            const double phase = std::real(potential[tick][i](c,c) * in.dt / 2.0);
            op[c][i][0] = cos(phase);
            op[c][i][1] = -sin(phase);
        }
    }
}
void definePotentialOperatorDiag(const inputs& in, const int& tick, const std::vector<fftw_complex*>& op,
    const std::function<hermitian_matrix(double, double)>& potential) {
    #pragma omp parallel for default(none) shared(in, potential, tick, op)
    for (int i = 0; i < in.space_grid; i++) {
        auto intermediate = potential(i*in.dx + in.initial_pos, (tick+0.5)*in.dt + in.initial_t);
        for (int c = 0; c < in.channels; c++) {
            const double phase = std::real(intermediate(c,c) * in.dt / 2.0);
            op[c][i][0] = cos(phase);
            op[c][i][1] = -sin(phase);
        }
    }
}
// creates array of diagonal potential operator arrays from potential at fractional tick and with coefficient
void definePotentialOperatorDiagFrac(const inputs& in, const double& tick, const std::vector<fftw_complex*>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential, const double& alpha) {
    #pragma omp parallel for default(none) shared(in, potential, tick, op, alpha)
    for (int i = 0; i < in.space_grid; i++) {
        for (int c = 0; c < in.channels; c++) {
            const double phase = alpha*std::real(potential[static_cast<int>(tick)][i](c,c) * in.dt);
            op[c][i][0] = cos(phase);
            op[c][i][1] = -sin(phase);
        }
    }
}
void definePotentialOperatorDiagFrac(const inputs& in, const double& tick, const std::vector<fftw_complex*>& op,
    const std::function<hermitian_matrix(double, double)>& potential, const double& alpha) {
    #pragma omp parallel for default(none) shared(in, potential, tick, op, alpha)
    for (int i = 0; i < in.space_grid; i++) {
        auto intermediate = potential(i*in.dx + in.initial_pos, tick*in.dt + in.initial_t);
        for (int c = 0; c < in.channels; c++) {
            const double phase = alpha*std::real(intermediate(c,c) * in.dt);
            op[c][i][0] = cos(phase);
            op[c][i][1] = -sin(phase);
        }
    }
}

// executes fft plans in parallel for all channels
void fftExecuteCC(const int gridpoints, const int channels, const std::vector<fftw_complex*>& psi,
    const std::vector<fftw_plan>& fft_plans) {
    #pragma omp parallel for default(none) shared(channels, gridpoints, psi, fft_plans)
    for (int c = 0; c < channels; c++) {
        for (int i = 0; i < gridpoints; i++) {
            const int sign = i % 2 == 0 ? 1 : -1;
            psi[c][i][0] *= sign;
            psi[c][i][1] *= sign;
        }
        fftw_execute(fft_plans[c]);
    }
}

// executes ifft plans in parallel for all channels
void ifftExecuteCC(const int gridpoints, const int channels, const std::vector<fftw_complex*>& psi,
    const std::vector<fftw_plan>& ifft_plans) {
    #pragma omp parallel for default(none) shared(channels, gridpoints, psi, ifft_plans)
    for (int c = 0; c < channels; c++) {
        fftw_execute(ifft_plans[c]);
        for (int i = 0; i < gridpoints; i++) {
            const int sign = i % 2 == 0 ? 1 : -1;
            psi[c][i][0] *= sign;
            psi[c][i][1] *= sign;
        }
    }
}

// applies diagonal part of potential operator
void applyDiagonalOperator(const int point, const int channels, const std::vector<fftw_complex*>& diag,
    const std::vector<fftw_complex*>& psi) {
    for (int c = 0; c < channels; c++) {
        const double re = psi[c][point][0];
        const double im = psi[c][point][1];
        psi[c][point][0] = re*diag[c][point][0] - im*diag[c][point][1];
        psi[c][point][1] = im*diag[c][point][0] + re*diag[c][point][1];
    }
}

// applies kinetic energy operator using applyKineticOperator from propagate.cpp in parallel
void applyKineticOperatorCC(const int gridpoints, const int channels, const std::vector<fftw_complex*>& psi,
    const fftw_complex* T) {
    #pragma omp parallel for default(none) shared(channels, gridpoints, psi, T)
    for (int c = 0; c < channels; c++) {
        applyKineticOperator(gridpoints, psi[c], T);
    }
}
// applies kinetic energy operator using applyKineticOperatorCoeff from propagate.cpp in parallel
void applyKineticOperatorCCCoeff(const int gridpoints, const int channels, const std::vector<fftw_complex*>& psi,
    const fftw_complex* T, const double& alpha, const double& dx, const double& dt) {
    #pragma omp parallel for default(none) shared(channels, gridpoints, psi, T, alpha, dx, dt)
    for (int c = 0; c < channels; c++) {
        applyKineticOperatorCoeff(gridpoints, psi[c], T, alpha, dx, dt);
    }
}

fftwResourcesCC fftwPrepTDCC(const inputs& in, const std::vector<fftw_complex*>& psi, const std::string& data) {
    const std::string wisdomfile = data + "/fftw_wisdom.dat";
    fftw_import_wisdom_from_filename(wisdomfile.c_str());
    std::vector<std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>> fft_ptrs;
    std::vector<std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>> ifft_ptrs;
    std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> V_d;
    for (int c = 0; c < in.channels; c++) {
        fft_ptrs.emplace_back(
            fftw_plan_dft_1d(in.space_grid, psi[c], psi[c], FFTW_FORWARD, FFTW_MEASURE),
            &fftw_destroy_plan
        );
        ifft_ptrs.emplace_back(
            fftw_plan_dft_1d(in.space_grid, psi[c], psi[c], FFTW_BACKWARD, FFTW_MEASURE),
            &fftw_destroy_plan
        );
        V_d.emplace_back(fftw_alloc_complex(in.space_grid), fftw_free);
    }
    fftw_export_wisdom_to_filename(wisdomfile.c_str());
    const auto T = fftw_alloc_complex(in.space_grid);
    defineKineticOperator(in, T, false);
    return fftwResourcesCC{
        .fft_ptrs = std::move(fft_ptrs),
        .ifft_ptrs = std::move(ifft_ptrs),
        .Tp = std::unique_ptr<fftw_complex, void(*)(void*)>(T, fftw_free),
        .V_d = std::move(V_d),
        .V_c = std::vector(in.space_grid, std::vector(in.channels, std::vector<coupling>(in.channels)))
    };
}

// normalizes fftw results in parallel for all channels
void fftwNormCC(const inputs& in, const double scale, const std::vector<fftw_complex*>& psi) {
    #pragma omp parallel for default(none) shared(in, scale, psi)
    for (int c = 0; c < in.channels; c++) {
        scale_fftw_complex(scale, psi[c], in.space_grid);
    }
}

// gets initial wavepacket handled with RAII, outputs a vector of unique pointers
std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> getWavepacketCC(const inputs& in, const std::string& data) {
    std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> psip;
    const auto psi = fftw_alloc_complex(in.space_grid);
    fftw_complex_array_from_file(data + "/" + in.input_psi_file + ".dat", psi, in.space_grid);
    psip.emplace_back(psi, fftw_free);
    for (int c = 1; c < in.channels; c++) {
        const auto psic = fftw_alloc_complex(in.space_grid);
        #pragma omp parallel for default(none) shared(in, psic)
        for (int i = 0; i < in.space_grid; i++) {
            psic[i][0] = 0.0;
            psic[i][1] = 0.0;
        }
        psip.emplace_back(psic, fftw_free);
    }
    return psip;
}

// writes output serially across channels (cannot be done in parallel)
void writeOutputCC(const std::vector<fftw_complex*>& psi, const int t, const inputs& in, std::vector<double>& buffer) {
    for (int c = 0; c < in.channels; c++) {
        writeOutput(psi[c], t, in.space_grid, in.nx_prints, in.nt_prints, buffer);
    }
}

// calculate and save PEY to file
void savePEY(const std::vector<fftw_complex*>& psi, const inputs& in) {
    std::vector<double> PEY(in.channels-1);
    double total_PEY = 0;
    for (int c = 1; c < in.channels; c++) {
        PEY[c-1] = norm(in.space_grid, in.dx, psi[c])/in.dE;
        total_PEY += PEY[c-1]*in.dE;
    }
    writeArray1D(in.min_energy, in.dE, in.channels-1, "data/PEY.dat", PEY);
    spacer(RESET);
    std::print("{}——> PEY saved! Total yield is {}{}{}\n", GRAY, YELLOW, total_PEY, RESET);
}
