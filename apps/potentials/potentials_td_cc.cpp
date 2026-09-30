//
// Created by Arian Dovald on 9/28/26.
//

#include "potentials_td_cc.h"
#include <numbers>

// harmonic oscillator coupled to decaying exponential by constant
double SHO(const double x) {
    return x*x/2;
}
double expDecay(const double x) {
    return 500*exp(-x);
}
double gaussian(const double x, const double delta, const double pos) {
    return pow(1/(std::numbers::pi * (delta*delta)), 1.0/4.0)*exp(-(x-pos)*(x-pos)/(2*(delta*delta)));
}
std::function<hermitian_matrix(double, double)> coupledSHO(const double coupling_strength) {
    return [coupling_strength](const double x, const double t) {
        hermitian_matrix potential(2);
        potential(0,0) = std::complex(SHO(x), 0.0);
        potential(1,1) = std::complex(expDecay(x), 0.0);
        potential(0,1) = std::complex(coupling_strength*gaussian(t, 0.5, 1), 0.0);
        return potential;
    };
}

// potential builder + map of options
std::function<hermitian_matrix(double, double)> buildPotentialTDCC(char* argv[]) {
    std::unordered_map<std::string,
    std::function<std::function<hermitian_matrix(double, double)>()>> const potentials = {
        {"coupledsho",    [argv] {
            return coupledSHO(std::stod(argv[6]));
        }}
    };
    return potentials.at(argv[3])();
}
