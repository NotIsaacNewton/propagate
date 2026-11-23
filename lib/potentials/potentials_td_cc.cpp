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

// gets all potentials
void potentials_from_file::getPotentials(const inputs& in) {
    for (int c = 0; c < in.bound_channels; c++) {
        readArray1D(in.pot_read_from_file + "_" + std::to_string(c) + ".dat", potentials[c]);
    }
}
// generates interpolators for all potentials
std::vector<spline_interp> potentials_from_file::interpPotentials(const inputs& in, const std::vector<double>& grid,
    const std::vector<std::vector<double>>& potentials) {
    std::vector<spline_interp> interps;
    for (int c = 0; c < in.bound_channels+in.residual_channels; c++) {
        spline_interp(grid, potentials[c]);
    }
    return interps;
}

// grabs all b-c dipole couplings into hermitian_matrix coup
void couplings_from_file::couplingsFromFileBtoC(const inputs& in) {
    std::vector<double> temp(2*in.residual_channels*in.continuum_channels);
    for (int c = 1; c < in.residual_channels; c++) {
        readArray1D(in.pot_read_from_file + "_0_" + std::to_string(c) + ".dat", temp);
    }
    for (int c = 1; c < in.residual_channels; c++) {
        coup(0,c) = {temp[2*(c-1)], temp[2*(c-1)+1]};
    }
}
// grabs a c-c dipole couplings ordered as (E_0->E_0, E_0->E_1, ..., E_0->E_M; ...; E_M->E_0, E_M->E_1, ... E_M->E_M)
void couplings_from_file::couplingsFromFileCtoC(const inputs& in) {
    std::vector temp(in.residual_channels, std::vector(in.residual_channels-1,
        std::vector<double>(2*in.continuum_channels*in.continuum_channels)));
    for (int c1 = 1; c1 < in.residual_channels; c1++) {
        for (int c2 = c1+1; c2 < in.residual_channels; c2++) {
            readArray1D(in.pot_read_from_file + "_" + std::to_string(c1) + "_" + std::to_string(c2) + ".dat",
                temp[c1][c2]);
        }
    }
    for (int c1 = 1; c1 < in.residual_channels; c1++) {
        for (int c2 = c1+1; c2 < in.residual_channels; c2++) {
            #pragma omp parallel for default(none) shared(in, temp, c1, c2)
            for (int i = 0; i < in.continuum_channels; i++) {
                for (int j = 0; j < in.continuum_channels; j++) {
                    coup(c1*in.continuum_channels + i, c2*in.residual_channels + j) =
                        {temp[c1][c2][2*i*in.continuum_channels + 2*j],
                        temp[c1][c2][2*i*in.continuum_channels + 2*j + 1]};
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
            }}
    };
    return potentials.at(in.potential_type)();
}
