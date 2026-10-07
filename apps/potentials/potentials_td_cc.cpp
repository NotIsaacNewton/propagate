//
// Created by Arian Dovald on 9/28/26.
//

#include "potentials_td_cc.h"
#include <numbers>
#include <string>
#include <cmath>
#include <algorithm>

// harmonic oscillator coupled to decaying exponential by constant
double SHO(const double x) {
    return x*x/2;
}
double expDecay(const double x, const int c) {
    return 500*exp(-x) + (c-1)*20;
}
double gaussian(const double x, const double delta, const double pos) {
    return pow(1/(std::numbers::pi * (delta*delta)), 1.0/4.0)*exp(-(x-pos)*(x-pos)/(2*(delta*delta)));
}
std::function<hermitian_matrix(double, double)> coupledSHO(const double coupling_strength, const int channels) {
    return [coupling_strength, channels](const double x, const double t) {
        hermitian_matrix potential(channels);
        potential(0,0) = std::complex(SHO(x), 0.0);
        for (int c = 1; c < channels; c++) {
            potential(c,c) = std::complex(expDecay(x, c), 0.0);
        }
        for (int c1 = 0; c1 < channels; c1++) {
            for (int c2 = c1+1; c2 < channels; c2++) {
                potential(c1,c2) = std::complex(coupling_strength*gaussian(t, 0.4, 2), 0.0);
            }
        }
        return potential;
    };
}

// potential builder + map of options
std::function<hermitian_matrix(double, double)> buildPotentialTDCC(const inputs& in) {
    std::unordered_map<std::string,
    std::function<std::function<hermitian_matrix(double, double)>()>> const potentials = {
        {"coupledsho",    [in] {
            return coupledSHO(in.strength_1, in.channels);
        }}
    };
    return potentials.at(in.potential_type)();
}
