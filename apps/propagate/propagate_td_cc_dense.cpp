//
// Created by Arian Dovald on 11/16/26.
//

#include "propagate_td_cc_dense.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>
#include <print>
#include <algorithm>
#include <variant>
#include <omp.h>
#include "../potentials/potentials_td_cc.h"
#include "../potentials/potentials_td.h"

// applies potential operator using dense hermitian matrix diagonalization ZHEEVR routine from LAPACK
void applyPotentialOperatorDense(const inputs& in, const int tick,
    const std::vector<std::vector<hermitian_matrix>>& potential, const std::vector<fftw_complex*>& psi) {
    // parallel over the spatial grid
    #pragma omp parallel for default(none) shared(in, tick, potential, psi)
    for (int i = 0; i < in.space_grid; i++) {
        const hermitian_matrix H = potential[tick][i];
        // allocate eigenvalue and eigenvector arrays
        std::vector<double> eigenvalues(in.channels);
        std::vector<lapack_complex_double> eigenvectors(in.channels * in.channels);
        std::vector<lapack_int> isuppz(2 * in.channels);
        // call ZHEEVR with error check
        diagonalizeHermitianMatrix(H, eigenvalues, eigenvectors, isuppz);
        // apply e^{-i * eigenvalue * dt/2} in eigenbasis
        // project psi onto eigenvectors
        std::vector<std::complex<double>> psi_vec(in.channels);
        for (int c = 0; c < in.channels; c++) {
            psi_vec[c] = {psi[c][i][0], psi[c][i][1]};
        }
        // coefficients
        std::vector<std::complex<double>> coeff(in.channels, {0.0, 0.0});
        for (int k = 0; k < in.channels; k++) {
            for (int c = 0; c < in.channels; c++) {
                const auto z = reinterpret_cast<const std::complex<double>*>(eigenvectors.data())[c + k*in.channels];
                coeff[k] += std::conj(z) * std::complex(psi_vec[c].real(), psi_vec[c].imag());
            }
        }
        // apply e^{-i * lambda[k] * dt/2}
        for (int k = 0; k < in.channels; k++) {
            const double phase = -eigenvalues[k] * in.dt / 2.0;
            coeff[k] *= std::complex(std::cos(phase), std::sin(phase));
        }
        // transform back: psi = Z * coeff
        for (int c = 0; c < in.channels; c++) {
            std::complex val = {0.0, 0.0};
            for (int k = 0; k < in.channels; k++) {
                const auto z = reinterpret_cast<const std::complex<double>*>(eigenvectors.data())[c + k*in.channels];
                val += z * coeff[k];
            }
            psi[c][i][0] = val.real();
            psi[c][i][1] = val.imag();
        }
    }
}
void applyPotentialOperatorDense(const inputs& in, const int tick,
    const std::function<hermitian_matrix(double,double)>& potential, const std::vector<fftw_complex*>& psi) {
    // parallel over the spatial grid
    #pragma omp parallel for default(none) shared(in, tick, potential, psi)
    for (int i = 0; i < in.space_grid; i++) {
        hermitian_matrix H = potential(i*in.dx + in.initial_pos, (tick+0.5)*in.dt + in.initial_t);
        for (int c1 = 0; c1 < in.channels; c1++) {
            for (int c2 = c1+1; c2 < in.channels; c2++) {
                H(c1, c2) *= (in.initial_pos + i*in.dx)*electricField(in.initial_t + (tick + 0.5)*in.dt,
                    in.strength_2);
            }
        }
        // allocate eigenvalue and eigenvector arrays
        std::vector<double> eigenvalues(in.channels);
        std::vector<lapack_complex_double> eigenvectors(in.channels * in.channels);
        std::vector<lapack_int> isuppz(2 * in.channels);
        // call ZHEEVR with error check
        diagonalizeHermitianMatrix(H, eigenvalues, eigenvectors, isuppz);
        // apply e^{-i * eigenvalue * dt/2} in eigenbasis
        // project psi onto eigenvectors
        std::vector<std::complex<double>> psi_vec(in.channels);
        for (int c = 0; c < in.channels; c++) {
            psi_vec[c] = {psi[c][i][0], psi[c][i][1]};
        }
        // coefficients = psi^* * psi
        std::vector<std::complex<double>> coeff(in.channels, {0.0, 0.0});
        for (int k = 0; k < in.channels; k++) {
            for (int c = 0; c < in.channels; c++) {
                const auto z = reinterpret_cast<const std::complex<double>*>(eigenvectors.data())[
                    c + k*in.channels];
                coeff[k] += std::conj(z) * std::complex(psi_vec[c].real(), psi_vec[c].imag());
            }
        }
        // apply e^{-i * lambda[k] * dt/2}
        for (int k = 0; k < in.channels; k++) {
            const double phase = -eigenvalues[k] * in.dt / 2.0;
            coeff[k] *= std::complex(std::cos(phase), std::sin(phase));
        }
        // transform back: psi = Z * coeff
        for (int c = 0; c < in.channels; c++) {
            std::complex val = {0.0, 0.0};
            for (int k = 0; k < in.channels; k++) {
                const auto z = reinterpret_cast<const std::complex<double>*>(eigenvectors.data())[
                    c + k*in.channels];
                val += z * coeff[k];
            }
            psi[c][i][0] = val.real();
            psi[c][i][1] = val.imag();
        }
    }
}

// propagates psi in potential from tick to tick + 1
void propTickDense(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, const double scale) {
    // apply potential operators
    std::visit([&](const auto& pot) {
        applyPotentialOperatorDense(in, tick, pot, psi);
    }, potential);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCC(in.space_grid, in.channels, psi, T);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // apply potential operators
    std::visit([&](const auto& pot) {
        applyPotentialOperatorDense(in, tick, pot, psi);
    }, potential);
}

// propagates wavefunction based on general values
void propagateTDCCDense(const inputs& in, const std::string& data) {
    // get initial wavepacket with RAII
    // psip owns the lifetime of psi's data
    auto psip = getWavepacketCC(in, data); // vector of pointers
    std::vector<fftw_complex*> psi(in.channels);
    for (int c = 0; c < in.channels; c++) {
        psi[c] = psip[c].get();
    }
    // prep fftw variables and plans
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    auto [fft_ptrs,
        ifft_ptrs,
        Tp,
        V_d_raw,
        V_c] = fftwPrepTDCC(in, psi, data);
    std::vector<fftw_complex*> V_d(in.channels);
    std::vector<fftw_plan> fft(in.channels), ifft(in.channels);
    for (int c = 0; c < in.channels; c++) {
        V_d[c] = V_d_raw[c].get();
        fft[c] = fft_ptrs[c].get();
        ifft[c] = ifft_ptrs[c].get();
    }
    fftw_complex* T = Tp.get();
    // get potential
    potCC potential = in.potential_type == "file"
        ? potCC(getPotentialCC(in, data))
        : potCC(buildPotentialTDCC(in));
    // scale for normalizing fft result
    const double scale = 1.0 / in.space_grid;
    // open output file and buffer
    auto [wf, buffer] = openWFOutputFile(in, data);
    // norm psi (often slightly off norm) and check norm
    normalizeCC(in.space_grid, in.dx, in.channels, psi); // normalize psi
    double mag = normCC(in.space_grid, in.dx, in.channels, psi); // recalculate norm
    std::print("The initial norm is {}\n", mag);
    // spacer
    spacer(RESET);
    // console output
    std::cout << "Propagating...\n";
    // propagation loop
    for (int t = 0; t < in.time_grid; t++) {
        // print completion % to console
        !((t+1) % (in.time_grid / 10)) ? progressBar(GREEN, 100*(t+1)/in.time_grid) : reset();
        // write lines in output buffer
        writeOutputCC(psi, t, in, buffer);
        // propagate for one tick
        propTickDense(t, in, psi, potential, T, fft, ifft, scale);
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
    mag = normCC(in.space_grid, in.dx, in.channels, psi);
    std::cout << RED << "Check norm:\n" << RESET;
    std::print("The final norm is {}\n", mag);
    // save PEY
    savePEY(psi, in);
}
