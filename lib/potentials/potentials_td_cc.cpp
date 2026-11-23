//
// Created by Arian Dovald on 9/28/26.
//

#include "potentials_td_cc.h"
#include <string>
#include <cmath>

// harmonic oscillator coupled to decaying exponential by constant
double SHO(const double x) {
    return x*x/2;
}
double expDecay(const double x, const int c, const int channels, const double range) {
    return 500*exp(-x) + (c-1)*(range/(channels-1));
}
std::function<hermitian_matrix(double, double)> coupledSHO(const inputs& in) {
    return [in](const double x, const double) {
        hermitian_matrix potential(in.channels);
        potential(0,0) = std::complex(SHO(x), 0.0);
        for (int c = 1; c < in.residual_channels*in.continuum_channels; c++) {
            potential(c,c) = std::complex(expDecay(x, c, in.continuum_channels,
                in.max_energy-in.min_energy), 0.0);
        }
        for (int c2 = 1; c2 < in.residual_channels*in.continuum_channels; c2++) {
            potential(0,c2) = {in.strength_1*exp(-(c2*in.dE-2)*(c2*in.dE-2)/(2*(0.5*0.5))), 0.0};
        }
        for (int c1 = 1; c1 < in.residual_channels*in.continuum_channels; c1++) {
            for (int c2 = c1+1; c2 < in.residual_channels*in.continuum_channels; c2++) {
                potential(c1,c2) = {in.strength_1*exp(-(c2*in.dE-2)*(c2*in.dE-2)/(2*(0.5*0.5))), 0.0};
            }
        }
        return potential;
    };
}

// this is called for the potentials during hybrid runs
std::function<std::vector<double>(double, int)> hybridPotentialFunctions(const int states, const double dE) {
    return [states, dE] (const double x, const int c) {
        std::vector<double> pots(states);
        pots[0] = -1/x; // a bound state
        pots[1] = -2/x + c*dE; // a residual state
        pots[2] = -3/x + c*dE; // another residual state
        return pots;
    };
}
std::function<hermitian_matrix(double, double)> hybridPotential(const inputs& in) {
    return [in](const double x, const double) {
        hermitian_matrix potential(in.bound_channels + in.residual_channels);
        const auto potentials =
            hybridPotentialFunctions(in.bound_channels + in.residual_channels, in.dE);
        potential(0,0) = std::complex(potentials(x, 0)[0], 0.0);
        for (int c = 1; c < in.residual_channels*in.continuum_channels; c++) {
            potential(c,c) = std::complex(expDecay(x, c, in.continuum_channels,
                in.max_energy-in.min_energy), 0.0);
        }
        return potential;
    };
}

// gets all potentials
void potentials_from_file::getPotentials(const inputs& in) {
    for (int c = 0; c < in.bound_channels+in.residual_channels; c++) {
        readArray1D(in.pot_read_from_file + "_" + std::to_string(c) + ".dat", potentials[c]);
    }
}
// generates interpolators for all potentials
std::vector<spline_interp> potentials_from_file::interpPotentials(const inputs& in, const std::vector<double>& grid,
    const std::vector<std::vector<double>>& potentials) {
    std::vector<spline_interp> interps;
    for (int c = 0; c < in.bound_channels+in.residual_channels; c++) {
        interps[c] = spline_interp(grid, potentials[c]);
    }
    return interps;
}

// grabs all b-c dipole couplings into hermitian_matrix coup
void couplings_from_file::couplingsFromFileBtoC(const inputs& in) {
    std::vector temp(in.residual_channels, std::vector<std::complex<double>>(in.continuum_channels));
    for (int c = 1; c < in.residual_channels; c++) {
        readComplexArray1D(in.pot_read_from_file + "_0_" + std::to_string(c) + ".dat",
            temp[c-1], in.continuum_channels);
    }
    for (int c = 1; c < in.residual_channels; c++) {
        for (int cont = 0; cont < in.continuum_channels; cont++) {
            coup(0,(c+1)*cont) = temp[c-1][cont];
        }
    }
}
// grabs a c-c dipole couplings ordered as (E_0->E_0, E_0->E_1, ..., E_0->E_M; ...; E_M->E_0, E_M->E_1, ... E_M->E_M)
void couplings_from_file::couplingsFromFileCtoC(const inputs& in) {
    std::vector temp(in.residual_channels, std::vector(in.residual_channels-1,
        std::vector<std::complex<double>>(in.continuum_channels*in.continuum_channels)));
    for (int c1 = 1; c1 < in.residual_channels; c1++) {
        for (int c2 = c1+1; c2 < in.residual_channels; c2++) {
            readComplexArray1D(in.pot_read_from_file + "_" + std::to_string(c1) + "_" + std::to_string(c2) + ".dat",
                temp[c1][c2-1], in.continuum_channels);
        }
    }
    for (int c1 = 1; c1 < in.residual_channels; c1++) {
        for (int c2 = c1+1; c2 < in.residual_channels; c2++) {
            #pragma omp parallel for default(none) shared(in, temp, c1, c2)
            for (int i = 0; i < in.continuum_channels; i++) {
                for (int j = 0; j < in.continuum_channels; j++) {
                    coup(in.bound_channels + c1*in.continuum_channels + i,
                        in.bound_channels + c2*in.continuum_channels + j) = temp[c1][c2-1][i*in.continuum_channels + j];
                }
            }
        }
    }
}

// potential builder + map of options
std::function<hermitian_matrix(double, double)> buildPotentialTDCC(const inputs& in) {
    std::unordered_map<std::string,
    std::function<std::function<hermitian_matrix(double, double)>()>> const potentials = {
        {"coupledsho",    [in] {
            return coupledSHO(in);
            }},
        {"file",            [in] {
            auto pot = std::make_shared<potential_from_file_tdcc>(in);
            return [pot, in](const double x, const double t) {
                return pot->operator()(x, t);
            };
            }},
        {"hybrid",            [in] {
            auto pot = std::make_shared<hybrid_from_file_tdcc>(in);
            return [pot, in](const double x, const double t) {
                return pot->operator()(x, t);
            };
            }}
    };
    return potentials.at(in.potential_type)();
}
