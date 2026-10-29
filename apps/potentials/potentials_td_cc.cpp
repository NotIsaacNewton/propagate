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
double expDecay(const double x, const int c, const int channels, const double range) {
    return 500*exp(-x) + (c-1)*(range/(channels-1));
}
double gaussian(const double x, const double t, const double delta, const double pos) {
    return pow(1/(std::numbers::pi * (delta*delta)), 1.0/4.0)*exp(-(x-(t-pos))*(x-(t-pos))/(2*(delta*delta)));
}
std::function<hermitian_matrix(double, double)> coupledSHO(const inputs& in) {
    return [in](const double x, const double t) {
        hermitian_matrix potential(in.channels);
        potential(0,0) = std::complex(SHO(x), 0.0);
        for (int c = 1; c < in.channels; c++) {
            potential(c,c) = std::complex(expDecay(x, c, in.channels,
                in.max_energy-in.min_energy), 0.0);
        }
        for (int c2 = 1; c2 < in.channels; c2++) {
            potential(0,c2) = std::complex(
                in.strength_1*gaussian(x, t, 0.4, 4)*cos(in.strength_2/137*x-in.strength_2*t), 0.0
                );
        }
        for (int c1 = 1; c1 < in.channels; c1++) {
            for (int c2 = c1+1; c2 < in.channels; c2++) {
                potential(c1,c2) = std::complex(
                    in.strength_1*gaussian(x, t, 0.4, 4)*cos(in.strength_2/137*x-in.strength_2*t), 0.0
                    );
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
            return coupledSHO(in);
        }}
    };
    return potentials.at(in.potential_type)();
}
