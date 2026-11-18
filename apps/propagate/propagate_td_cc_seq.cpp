//
// Created by Arian Dovald on 11/16/26.
//

#include "propagate_td_cc_seq.h"
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

// creates array of coupling operator arrays from data at tick and outputs to op
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<std::vector<coupling>>>& op,
    const std::function<hermitian_matrix(double,double)>& potential) {
    #pragma omp parallel for default(none) shared(in, potential, tick, op)
    for (int i = 0; i < in.space_grid; i++) {
        auto temp = potential(i*in.dx + in.initial_pos, (tick+0.5)*in.dt + in.initial_t);
        for (int c1 = 0; c1 < in.channels; c1++) {
            for (int c2 = c1+1; c2 < in.channels; c2++) {
                const std::complex<double> coupling = (in.initial_pos + i*in.dx)*electricField(
                    in.initial_t + (tick + 0.5)*in.dt, in.strength_2) * temp(c1,c2);
                if (const double coupling_strength = std::abs(coupling); coupling_strength < 1e-15) {
                    op[i][c1][c2].cos_factor = 1.0;
                    op[i][c1][c2].sin_factor = 0.0;
                    op[i][c1][c2].exp_phase = std::complex(0.0, 0.0);
                } else {
                    const double phase = coupling_strength * in.dt / 2.0;
                    op[i][c1][c2].cos_factor = cos(phase);
                    op[i][c1][c2].sin_factor = sin(phase);
                    op[i][c1][c2].exp_phase = std::complex(coupling.imag()/coupling_strength,
                        -coupling.real()/coupling_strength);
                }
            }
        }
    }
}
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<std::vector<coupling>>>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential) {
    #pragma omp parallel for default(none) shared(in, potential, tick, op)
    for (int i = 0; i < in.space_grid; i++) {
        auto intermediate = potential[tick][i];
        for (int c1 = 0; c1 < in.channels; c1++) {
            for (int c2 = c1+1; c2 < in.channels; c2++) {
                const std::complex<double> coupling = intermediate(c1,c2);
                if (const double coupling_strength = std::abs(coupling); coupling_strength < 1e-15) {
                    op[i][c1][c2].cos_factor = 1.0;
                    op[i][c1][c2].sin_factor = 0.0;
                    op[i][c1][c2].exp_phase = std::complex(0.0, 0.0);
                } else {
                    const double phase = coupling_strength * in.dt / 2.0;
                    op[i][c1][c2].cos_factor = cos(phase);
                    op[i][c1][c2].sin_factor = sin(phase);
                    op[i][c1][c2].exp_phase = std::complex(coupling.imag()/coupling_strength,
                        -coupling.real()/coupling_strength);
                }
            }
        }
    }
}

// applies coupling part of potential operator for all channels at a specific gridpoint
void applyCouplingOperator(const int point, const int channels,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi) {
    for (int c1 = 0; c1 < channels; c1++) {
        for (int c2 = c1+1; c2 < channels; c2++) {
            const auto psi1_old = fftw_complex_to_std_complex(psi[c1][point]);
            const auto psi2_old = fftw_complex_to_std_complex(psi[c2][point]);
            std_complex_to_fftw_complex(psi1_old*coup[point][c1][c2].cos_factor +
                psi2_old*coup[point][c1][c2].exp_phase*coup[point][c1][c2].sin_factor, psi[c1][point]);
            std_complex_to_fftw_complex(psi2_old*coup[point][c1][c2].cos_factor -
                psi1_old*std::conj(coup[point][c1][c2].exp_phase)*coup[point][c1][c2].sin_factor, psi[c2][point]);
        }
    }
}
void applyCouplingOperatorReverse(const int point, const int channels,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi) {
    for (int c1 = channels-2; c1 >= 0; c1--) {
        for (int c2 = channels-1; c2 > c1; c2--) {
            const auto psi1_old = fftw_complex_to_std_complex(psi[c1][point]);
            const auto psi2_old = fftw_complex_to_std_complex(psi[c2][point]);
            std_complex_to_fftw_complex(psi1_old*coup[point][c1][c2].cos_factor +
                psi2_old*coup[point][c1][c2].exp_phase*coup[point][c1][c2].sin_factor, psi[c1][point]);
            std_complex_to_fftw_complex(psi2_old*coup[point][c1][c2].cos_factor -
                psi1_old*std::conj(coup[point][c1][c2].exp_phase)*coup[point][c1][c2].sin_factor, psi[c2][point]);
        }
    }
}

// applies potential operator in parallel threads across spatial gridpoints
void applyPotentialOperatorCC(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi) {
    #pragma omp parallel for default(none) shared(in, coup, psi, diag)
    for (int i = 0; i < in.space_grid; i++) {
        applyCouplingOperator(i, in.channels, coup, psi);
        applyDiagonalOperator(i, in.channels, diag, psi);
    }
}
void applyPotentialOperatorCCReverse(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi) {
    #pragma omp parallel for default(none) shared(in, coup, psi, diag)
    for (int i = 0; i < in.space_grid; i++) {
        applyDiagonalOperator(i, in.channels, diag, psi);
        applyCouplingOperatorReverse(i, in.channels, coup, psi);
    }
}

// propagates psi in potential from tick to tick + 1
void propTickTDCCSeq(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const std::vector<fftw_complex*>& V_d,
    std::vector<std::vector<std::vector<coupling>>>& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, const double scale) {
    // define diagonal and coupling potential operators
    std::visit([&](const auto& pot) {
        definePotentialOperatorDiag(in, tick, V_d, pot);
        defineCouplingOperator(in, tick, V_c, pot);
    }, potential);
    // apply potential operators
    applyPotentialOperatorCC(in, V_d, V_c, psi);
    // execute fft
    fftExecuteCC(in.space_grid, in.channels, psi, fft);
    // apply kinetic energy operator
    applyKineticOperatorCC(in.space_grid, in.channels, psi, T);
    // execute ifft
    ifftExecuteCC(in.space_grid, in.channels, psi, ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    fftwNormCC(in, scale, psi);
    // apply potential operators
    applyPotentialOperatorCCReverse(in, V_d, V_c, psi);
}

// propagates wavefunction based on general values
void propagateTDCCSeq(const inputs& in, const std::string& data) {
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
        propTickTDCCSeq(t, in, psi, potential, V_d, V_c, T, fft, ifft, scale);
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
