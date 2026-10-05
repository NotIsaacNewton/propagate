//
// Created by Arian Dovald on 8/20/26.
//

#include "potentials_td.h"
#include "potentials.h"
#include <cmath>
#include <string>

// well + wave
std::function<double(double, double)> electricBarrier(double start_pos, double end_pos,
    double strength_1, double strength_2) {
    return [start_pos, end_pos, strength_1, strength_2](const double x, const double t) {
        return x <= end_pos && x >= start_pos ?
        strength_1 + std::abs(strength_1)*cos(strength_2*x-137*strength_2*t)*exp(-(x-t+8)*(x-t+8)/8) :
        std::abs(strength_1)*cos(strength_2*x-137*strength_2*t)*exp(-(x-t+8)*(x-t+8)/8);
    };
}

// harmonic oscillator + wave
std::function<double(double, double)> electricSho(const double strength, const double k) {
    return [strength, k](const double x, const double t) {
        return strength*x*x/2 + std::abs(strength)*50*cos(k*x-137*k*t)*exp(-(x-t+4)*(x-t+4)/1);
    };
}

// potential builder + map of options
std::function<double(double, double)> buildPotentialTD(const inputs& in) {
    std::unordered_map<std::string, std::function<std::function<double(double, double)>()>> const potentials = {
        {"electricwell",    [in] {
            return electricBarrier(in.pot_position_1,in.pot_position_2,
                in.strength_1, in.strength_2);
            }},
        {"electricsho",    [in] {
            return electricSho(in.strength_1,in.strength_2);
            }}
    };
    return potentials.at(in.potential_type)();
}
