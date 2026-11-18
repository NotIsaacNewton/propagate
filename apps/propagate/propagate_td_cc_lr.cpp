//
// Created by Arian Dovald on 11/16/26.
//

#include "propagate_td_cc_lr.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <print>
#include <algorithm>
#include <variant>
#include <omp.h>
#include "../potentials/potentials_td_cc.h"
#include "../potentials/potentials_td.h"

// prep low-rank approximation method
coupling_lr prepCouplingOperatorLR(const inputs& in,
    const std::function<hermitian_matrix(double,double)>& potential) {
    hermitian_matrix coupling_matrix(in.channels);
    // grab dipole coupling matrix
    #pragma omp parallel for default(none) shared(in, potential, coupling_matrix)
    for (int c1 = 0; c1 < in.channels; c1++) {
        for (int c2 = c1+1; c2 < in.channels; c2++) {
            coupling_matrix(c1,c2) = potential(0.0,0.0)(c1,c2);
        }
    }
    // allocate eigenvalue and eigenvector arrays
    std::vector<double> eigenvalues(in.channels);
    std::vector<lapack_complex_double> eigenvectors(in.channels * in.channels);
    std::vector<lapack_int> isuppz(2 * in.channels);
    // diagonalize dipole coupling matrix
    diagonalizeHermitianMatrix(coupling_matrix, eigenvalues, eigenvectors, isuppz);
    const diag_matrix coup(in.channels, eigenvalues);
    return { .diag_coup = coup, .eigenvectors = eigenvectors };
}
coupling_lr prepCouplingOperatorLR(const inputs& in,
    const std::vector<std::vector<hermitian_matrix>>& potential) {
    hermitian_matrix coupling_matrix(in.channels);
    // grab dipole coupling matrix
    #pragma omp parallel for default(none) shared(in, potential, coupling_matrix)
    for (int c1 = 0; c1 < in.channels; c1++) {
        for (int c2 = c1+1; c2 < in.channels; c2++) {
            coupling_matrix(c1,c2) = potential[0][0](c1,c2);
        }
    }
    // allocate eigenvalue and eigenvector arrays
    std::vector<double> eigenvalues(in.channels);
    std::vector<lapack_complex_double> eigenvectors(in.channels * in.channels);
    std::vector<lapack_int> isuppz(2 * in.channels);
    // diagonalize dipole coupling matrix
    diagonalizeHermitianMatrix(coupling_matrix, eigenvalues, eigenvectors, isuppz);
    const diag_matrix coup(in.channels, eigenvalues);
    return { .diag_coup = coup, .eigenvectors = eigenvectors };
}

// apply low-rank coupling operator
void applyCouplingOperatorLR(const inputs& in, const coupling_lr& coup,
    const std::vector<fftw_complex*>& psi, const int point, const int tick) {
    // project psi onto eigenvectors
    std::vector<std::complex<double>> psi_vec(in.channels);
    for (int c = 0; c < in.channels; c++) { psi_vec[c] = {psi[c][point][0], psi[c][point][1]}; }
    // stores coefficients
    std::vector<std::complex<double>> coeff(in.channels);
    // apply e^{-i * x E(t) * lambda[k] * dt/2}
    // phases for transform
    std::vector<double> phases(in.channels);
    for (int k = 0; k < in.channels; k++) {
        phases[k] = - (in.initial_pos + point*in.dx)
        * electricField(in.initial_t + (tick+0.5)*in.dt, in.strength_2)
        * coup.diag_coup(k) * in.dt / 2.0;
    }
    // apply transformation in eigenbasis
    eigenbasisTransform(in.channels, coup.eigenvectors, phases, psi_vec, coeff);
    // write to psi
    for (int c = 0; c < in.channels; c++) {
        psi[c][point][0] = psi_vec[c].real();
        psi[c][point][1] = psi_vec[c].imag();
    }
}

// applies potential operator in parallel threads across spatial gridpoints for low-rank approximation
void applyPotentialOperatorLR(const inputs& in, const int tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi) {
    #pragma omp parallel for default(none) shared(in, coup, psi, diag, tick)
    for (int i = 0; i < in.space_grid; i++) {
        applyCouplingOperatorLR(in, coup, psi, i, tick);
        applyDiagonalOperator(i, in.channels, diag, psi);
    }
}
void applyPotentialOperatorLRReverse(const inputs& in, const int tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi) {
    #pragma omp parallel for default(none) shared(in, coup, psi, diag, tick)
    for (int i = 0; i < in.space_grid; i++) {
        applyDiagonalOperator(i, in.channels, diag, psi);
        applyCouplingOperatorLR(in, coup, psi, i, tick);
    }
}

// propagates psi in potential from tick to tick + 1
void propTickLR(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const std::vector<fftw_complex*>& V_d,
    const coupling_lr& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, const double scale) {
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiag(in, tick, V_d, pot);
    }, potential);
    // apply potential operators
    applyPotentialOperatorLR(in, tick, V_d, V_c, psi);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCC(in.space_grid, in.channels, psi, T);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // apply potential operators
    applyPotentialOperatorLRReverse(in, tick, V_d, V_c, psi);
}

// propagates wavefunction based on general values
void propagateTDCCLR(const inputs& in, const std::string& data) {
    // get initial wavepacket with RAII
    // psip owns the lifetime of psi's data
    auto psip = getWavepacketCC(in, data); // vector of pointers
    std::vector<fftw_complex*> psi(in.channels);
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
    // write to psi *after* FFTW_MEASURE is done (fftw planning can overwrite arrays)
    for (int c = 0; c < in.channels; c++) {
        psi[c] = psip[c].get();
    }
    // get potential
    potCC potential = in.potential_type == "file"
        ? potCC(getPotentialCC(in, data))
        : potCC(buildPotentialTDCC(in));
    coupling_lr coup = in.potential_type == "file"
        ? prepCouplingOperatorLR(in, getPotentialCC(in, data))
        : prepCouplingOperatorLR(in, buildPotentialTDCC(in));
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
        propTickLR(t, in, psi, potential, V_d, coup, T, fft, ifft, scale);
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
