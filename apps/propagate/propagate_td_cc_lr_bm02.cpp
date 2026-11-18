//
// Created by Arian Dovald on 11/17/26.
//

#include "propagate_td_cc_lr_bm02.h"
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

// apply low-rank coupling operator with modified coefficient alpha and with fractional (double-valued) tick
void applyCouplingOperatorLRCoeff(const inputs& in, const coupling_lr& coup,
    const std::vector<fftw_complex*>& psi, const int point, const double tick, const double alpha) {
    // project psi onto eigenvectors
    std::vector<std::complex<double>> psi_vec(in.channels);
    for (int c = 0; c < in.channels; c++) { psi_vec[c] = {psi[c][point][0], psi[c][point][1]}; }
    // stores coefficients
    std::vector<std::complex<double>> coeff(in.channels);
    // apply e^{-i * alpha * x E(t) * lambda[k] * dt}
    // phases for transform
    std::vector<double> phases(in.channels);
    for (int k = 0; k < in.channels; k++) {
        phases[k] = - alpha * (in.initial_pos + point*in.dx)
        * electricField(in.initial_t + tick*in.dt, in.strength_2)
        * coup.diag_coup(k) * in.dt;
    }
    // apply transformation in eigenbasis
    eigenbasisTransform(in.channels, coup.eigenvectors, phases, psi_vec, coeff);
    // write to psi
    for (int c = 0; c < in.channels; c++) {
        psi[c][point][0] = psi_vec[c].real();
        psi[c][point][1] = psi_vec[c].imag();
    }
}

// applies potential operator for low-rank approximation with coefficient alpha and fractional tick
void applyPotentialOperatorLRCoeff(const inputs& in, const double tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi, const double alpha) {
    #pragma omp parallel for default(none) shared(in, coup, psi, diag, tick, alpha)
    for (int i = 0; i < in.space_grid; i++) {
        applyCouplingOperatorLRCoeff(in, coup, psi, i, tick, alpha);
        applyDiagonalOperator(i, in.channels, diag, psi);
    }
}
void applyPotentialOperatorLRCoeffReverse(const inputs& in, const double tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi, const double alpha) {
    #pragma omp parallel for default(none) shared(in, coup, psi, diag, tick, alpha)
    for (int i = 0; i < in.space_grid; i++) {
        applyDiagonalOperator(i, in.channels, diag, psi);
        applyCouplingOperatorLRCoeff(in, coup, psi, i, tick, alpha);
    }
}

// propagates psi in potential from tick to tick + 1
void propTickLRBM02(const double& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const std::vector<fftw_complex*>& V_d,
    const coupling_lr& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, const double scale,
    const std::vector<double>& c) {
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiagFrac(in, tick, V_d, pot, c[0]);
    }, potential);
    // apply potential operators
    applyPotentialOperatorLRCoeff(in, tick, V_d, V_c, psi, c[0]);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCCCoeff(in.space_grid, in.channels, psi, T, c[1], in.dx, in.dt);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiagFrac(in, tick+c[1], V_d, pot, c[2]);
    }, potential);
    // apply potential operators
    applyPotentialOperatorLRCoeff(in, tick+c[1], V_d, V_c, psi, c[2]);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCCCoeff(in.space_grid, in.channels, psi, T, c[3], in.dx, in.dt);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiagFrac(in, tick+c[1]+c[3], V_d, pot, c[4]);
    }, potential);
    // apply potential operators
    applyPotentialOperatorLRCoeff(in, tick+c[1]+c[3], V_d, V_c, psi, c[4]);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCCCoeff(in.space_grid, in.channels, psi, T, c[5], in.dx, in.dt);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiagFrac(in, tick+c[1]+c[3]+c[5], V_d, pot, c[6]);
    }, potential);
    #pragma omp parallel for default(none) shared(in, V_c, psi, V_d, tick, c)
    for (int i = 0; i < in.space_grid; i++) {
        applyCouplingOperatorLRCoeff(in, V_c, psi, i, tick+c[1]+c[3]+c[5], c[6]/2.0);
        applyDiagonalOperator(i, in.channels, V_d, psi);
        applyCouplingOperatorLRCoeff(in, V_c, psi, i, tick+c[1]+c[3]+c[5], c[6]/2.0);
    }
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCCCoeff(in.space_grid, in.channels, psi, T, c[5], in.dx, in.dt);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiagFrac(in, tick+c[1]+c[3]+2*c[5], V_d, pot, c[4]);
    }, potential);
    // apply potential operators
    applyPotentialOperatorLRCoeffReverse(in, tick+c[1]+c[3]+2*c[5], V_d, V_c, psi, c[4]);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCCCoeff(in.space_grid, in.channels, psi, T, c[3], in.dx, in.dt);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiagFrac(in, tick+c[1]+2*c[3]+2*c[5], V_d, pot, c[2]);
    }, potential);
    // apply potential operators
    applyPotentialOperatorLRCoeffReverse(in, tick+c[1]+2*c[3]+2*c[5], V_d, V_c, psi, c[2]);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCCCoeff(in.space_grid, in.channels, psi, T, c[1], in.dx, in.dt);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // define diagonal potential operator
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiagFrac(in, tick+1, V_d, pot, c[0]);
    }, potential);
    // apply potential operators
    applyPotentialOperatorLRCoeffReverse(in, tick+1, V_d, V_c, psi, c[0]);
}

// propagates wavefunction based on general values
void propagateTDCCLRBM02(const inputs& in, const std::string& data) {
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
    std::vector bm02_coefficients = {0.0829844064174052, 0.2452989571842710, 0.396309801498368,
        0.604872665711080, -0.0390563049223486, -0.350171622895351, 0.11952419401315084};
    progressBar(GREEN, 0);
    for (int t = 0; t < in.time_grid; t++) {
        // print completion % to console
        !((t+1) % (in.time_grid / 10)) ? progressBar(GREEN, 100*(t+1)/in.time_grid) : reset();
        // write lines in output buffer
        writeOutputCC(psi, t, in, buffer);
        // propagate for one tick
        propTickLRBM02(t, in, psi, potential, V_d, coup, T, fft, ifft, scale, bm02_coefficients);
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
