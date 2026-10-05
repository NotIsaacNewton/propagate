//
// Created by Arian Dovald on 8/20/26.
//

#ifndef PROPAGATE_POTENTIALS_TD_H
#define PROPAGATE_POTENTIALS_TD_H

#include <functional>
#include "file_tools.h"

// well + wave
std::function<double(double, double)> electricBarrier(double start_pos, double end_pos,
    double strength_1, double strength_2);
// harmonic oscillator + wave
std::function<double(double, double)> electricSho(double strength, double k);

// potential builder + map of options
std::function<double(double, double)> buildPotentialTD(const inputs& in);

#endif //PROPAGATE_POTENTIALS_TD_H
