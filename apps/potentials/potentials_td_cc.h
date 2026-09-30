//
// Created by Arian Dovald on 9/28/26.
//

#ifndef PROPAGATE_POTENTIALS_TD_CC_H
#define PROPAGATE_POTENTIALS_TD_CC_H

#include <functional>
#include "matrix_tools.h"


// TODO:
//  some test potentials for continuum discretization
//  helper functions for reading in from a variety of file types

// harmonic oscillator coupled to decaying exponential by constant
double SHO(double x);
double expDecay(double x);
double gaussian(double x, double delta, double pos);
std::function<hermitian_matrix(double, double)> coupledSHO(double coupling_strength);

// potential builder + map of options
std::function<hermitian_matrix(double, double)> buildPotentialTDCC(char* argv[]);

#endif //PROPAGATE_POTENTIALS_TD_CC_H
