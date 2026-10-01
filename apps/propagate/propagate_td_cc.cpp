//
// Created by Arian Dovald on 9/2/26.
//

#include "propagate_td_cc.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"
#include <iostream>

// gets potential and returns arrays
std::vector<std::vector<hermitian_matrix> > getPotentialCC(const inputs &in, const std::string &data) {
    // prepare potential arrays
    std::vector potential(in.time_grid,
    std::vector(in.space_grid, hermitian_matrix(in.channels))); // full potential array
    std::vector<std::vector<double>> temp; // array containing raw values
    std::print("Reading potential from {}/potential.dat...\n", data);
    readArray2D(data + "/potential.dat", temp,
                in.space_grid*(in.channels*in.channels+in.channels),
                in.time_grid); // reads in potential file
    #pragma omp parallel for
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
    std::print("\nPotential read!\n", data);
    return potential;
}

// creates array of diagonal potential operator arrays from potential at tick and outputs to op
void definePotentialOperatorDiag(const inputs& in, const int& tick, const std::vector<fftw_complex*>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential) {
    for (int i = 0; i < in.space_grid; i++) {
        for (int c = 0; c < in.channels; c++) {
            const double phase = std::real(potential[tick][i](c,c) * in.dt / 2.0);
            op[c][i][0] = cos(phase);
            op[c][i][1] = -sin(phase);
        }
    }
}

// creates array of coupling operator arrays from data at tick and outputs to op
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<std::vector<coupling>>>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential) {
    for (int i = 0; i < in.space_grid; i++) {
        for (int c1 = 0; c1 < in.channels; c1++) {
            for (int c2 = c1+1; c2 < in.channels; c2++) {
                const std::complex<double> coupling = potential[tick][i](c1,c2);
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
    #pragma omp parallel for
    for (int i = 0; i < in.space_grid; i++) {
        applyCouplingOperator(i, in.channels, coup, psi);
        applyDiagonalOperator(i, in.channels, diag, psi);
    }
}
void applyPotentialOperatorCCReverse(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi) {
    #pragma omp parallel for
    for (int i = 0; i < in.space_grid; i++) {
        applyDiagonalOperator(i, in.channels, diag, psi);
        applyCouplingOperator(i, in.channels, coup, psi);
    }
}

// executes fft plans in parallel for all channels
void fftExecuteCC(const int gridpoints, const int channels, const std::vector<fftw_complex*>& psi,
    const std::vector<fftw_plan>& fft_plans) {
    #pragma omp parallel for
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
    #pragma omp parallel for
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
    #pragma omp parallel for
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
    #pragma omp parallel for
    for (int c = 0; c < in.channels; c++) {
        scale_fftw_complex(scale, psi[c], in.space_grid);
    }
}

// propagates psi in potential from tick to tick + 1
void propTickTDCC(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const std::vector<std::vector<hermitian_matrix>>& potential, const std::vector<fftw_complex*>& V_d,
    std::vector<std::vector<std::vector<coupling>>>& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, const double scale) {
    // define diagonal and coupling potential operators
    definePotentialOperatorDiag(in, tick, V_d, potential);
    defineCouplingOperator(in, tick, V_c, potential);
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

// gets initial wavepacket handled with RAII, outputing a vector of unique pointers
std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> getWavepacketCC(const inputs& in, const std::string& data) {
    std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> psip;
    for (int c = 0; c < in.channels; c++) {
        const auto psi = fftw_alloc_complex(in.space_grid);
        fftw_complex_array_from_file(data + "/psi_initial_" + std::to_string(c) + ".dat", psi, in.space_grid);
        psip.emplace_back(psi, fftw_free);
    }
    return psip;
}

// writes output serially across channels (cannot be done in parallel)
void writeOutputCC(const std::vector<fftw_complex*>& psi, const int t, const inputs& in, std::vector<double>& buffer) {
    for (int c = 0; c < in.channels; c++) {
        writeOutput(psi[c], t, in.space_grid, in.nx_prints, in.nt_prints, buffer);
    }
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
    auto potential = getPotentialCC(in, data);
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
    for (int t = 0; t < in.time_grid; t++) {
        // print completion % to console
        !((t+1) % (in.time_grid / 10)) ? progressBar(GREEN, 100*(t+1)/in.time_grid) : reset();
        // write lines in output file
        writeOutputCC(psi, t, in, buffer);
        // propagate for one tick
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
}
