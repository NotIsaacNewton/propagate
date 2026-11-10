//
// Created by Arian Dovald on 9/2/26.
//

#include "propagate_td_cc.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"
#include "lapacke.h"
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

// creates array of coupling operator arrays from data at tick and outputs to op
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<std::vector<coupling>>>& op,
    const std::function<hermitian_matrix(double,double)>& potential) {
    #pragma omp parallel for default(none) shared(in, potential, tick, op)
    for (int i = 0; i < in.space_grid; i++) {
        auto intermediate = potential(i*in.dx + in.initial_pos, (tick+0.5)*in.dt + in.initial_t);
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

// applies potential operator using dense hermitian matrix diagonalization ZHEEVR routine from LAPACK
void applyPotentialOperatorDense(const inputs& in, const int tick,
    const std::vector<std::vector<hermitian_matrix>>& potential, const std::vector<fftw_complex*>& psi) {
    // parallel over the spatial grid
    #pragma omp parallel for default(none) shared(in, tick, potential, psi)
    for (int i = 0; i < in.space_grid; i++) {
        const hermitian_matrix H = potential[tick][i];
        // unpack hermitian_matrix into column-major dense array for LAPACKE
        // NOTE: lapack_complex_double* is ABI-compatible with std::complex<double>
        std::vector<lapack_complex_double> A(in.channels * in.channels, {0.0, 0.0});
        for (int col = 0; col < in.channels; col++) {
            for (int row = 0; row <= col; row++) {
                // index = row + col*channels
                const auto val = H(row, col);
                A[row + col * in.channels] = {val.real(), val.imag()};
            }
        }

        // allocate eigenvalue and eigenvector arrays
        std::vector<double> eigenvalues(in.channels);
        std::vector<lapack_complex_double> eigenvectors(in.channels * in.channels);
        std::vector<lapack_int> isuppz(2 * in.channels);
        lapack_int m = 0;

        // call ZHEEVR
        lapack_int info = LAPACKE_zheevr(
            LAPACK_COL_MAJOR,           // storage order
            'V',                        // compute eigenvectors
            'A',                        // all eigenvalues
            'U',                        // upper triangle stored
            in.channels,                // matrix size
            A.data(),                   // matrix (overwritten)
            in.channels,                // leading dimension
            0.0, 0.0,                // vl, vu (unused for range='A')
            0, 0,                     // il, iu (unused for range='A')
            LAPACKE_dlamch('S'),  // abstol: safe minimum for reliability
            &m,                         // number of eigenvalues found
            eigenvalues.data(),
            eigenvectors.data(),
            in.channels,                // ldz
            isuppz.data()
        );
        // error check
        if (info != 0) {
            std::print("ZHEEVR failed at spatial point {} with info={}\n", i, info);
            continue;
        }

        // apply e^{-i * eigenvalue * dt/2} in eigenbasis
        // project psi onto eigenvectors
        std::vector<lapack_complex_double> psi_vec(in.channels);
        for (int c = 0; c < in.channels; c++) {
            psi_vec[c] = {psi[c][i][0], psi[c][i][1]};
        }
        // coefficients = psi^* * psi
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
        const hermitian_matrix H = potential(i*in.dx + in.initial_pos, (tick+0.5)*in.dt + in.initial_t);
        // unpack hermitian_matrix into column-major dense array for LAPACKE
        // NOTE: lapack_complex_double* is ABI-compatible with std::complex<double>
        std::vector<lapack_complex_double> A(in.channels * in.channels, {0.0, 0.0});
        for (int col = 0; col < in.channels; col++) {
            for (int row = 0; row <= col; row++) {
                // index = row + col*channels
                const auto val = H(row, col);
                A[row + col * in.channels] = {val.real(), val.imag()};
                // lower triangle (conjugate) — needed if uplo='U' isn't enough
                // A[col + row * N] = {val.real(), -val.imag()};
            }
        }

        // allocate eigenvalue and eigenvector arrays
        std::vector<double> eigenvalues(in.channels);
        std::vector<lapack_complex_double> eigenvectors(in.channels * in.channels);
        std::vector<lapack_int> isuppz(2 * in.channels);
        lapack_int m = 0;

        // call ZHEEVR
        lapack_int info = LAPACKE_zheevr(
            LAPACK_COL_MAJOR,           // storage order
            'V',                        // compute eigenvectors
            'A',                        // all eigenvalues
            'U',                        // upper triangle stored
            in.channels,                // matrix size
            A.data(),                   // matrix (overwritten)
            in.channels,                // leading dimension
            0.0, 0.0,                // vl, vu (unused for range='A')
            0, 0,                     // il, iu (unused for range='A')
            LAPACKE_dlamch('S'),  // abstol: safe minimum for reliability
            &m,                         // number of eigenvalues found
            eigenvalues.data(),
            eigenvectors.data(),
            in.channels,                // ldz
            isuppz.data()
        );
        // error check
        if (info != 0) {
            std::print("ZHEEVR failed at spatial point {} with info={}\n", i, info);
            continue;
        }

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

// applies kinetic energy operator using applyKineticOperator from propagate.cpp in parallel
void applyKineticOperatorCC(const int gridpoints, const int channels, const std::vector<fftw_complex*>& psi,
    const fftw_complex* T) {
    #pragma omp parallel for default(none) shared(channels, gridpoints, psi, T)
    for (int c = 0; c < channels; c++) {
        applyKineticOperator(gridpoints, psi[c], T);
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

// propagates psi in potential from tick to tick + 1
void propTickTDCC(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
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
    /* for (int c = 0; c < in.channels; c++) {
        const auto psi = fftw_alloc_complex(in.space_grid);
        fftw_complex_array_from_file(data + "/" + in.input_psi_file + "_" + std::to_string(c) + ".dat",
            psi, in.space_grid);
        psip.emplace_back(psi, fftw_free);
    } */
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
    std::print("\nPEY saved! Total yield is {}\n", total_PEY);
}

// propagates wavefunction based on general values
void propagateTDCC(const inputs& in, const std::string& data) {
    // get initial wavepacket with RAII
    // NOTE: psip owns the lifetime of psi's data
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
    std::cout << RED << "Check norm:\n" << RESET;
    std::print("The initial norm is {}\n", mag);
    // spacer
    spacer(RESET);
    // console output
    std::cout << "Propagating...\n";
    // propagation loop
    bool dense = in.run_type=="dense";
    for (int t = 0; t < in.time_grid; t++) {
        // print completion % to console
        !((t+1) % (in.time_grid / 10)) ? progressBar(GREEN, 100*(t+1)/in.time_grid) : reset();
        // write lines in output buffer
        writeOutputCC(psi, t, in, buffer);
        // propagate for one tick
        dense?  propTickDense(t, in, psi, potential, T, fft, ifft, scale) :
        propTickTDCC(t, in, psi, potential, V_d, V_c, T, fft, ifft, scale);
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
