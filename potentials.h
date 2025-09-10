//
// Created by Arian Dovald on 9/8/25.
//

#ifndef PROPAGATE_POTENTIALS_H
#define PROPAGATE_POTENTIALS_H

#include <functional>
#include "fftw3.h"

// creates potential operator array from potential function and outputs to op
[[maybe_unused]] inline void definePotentialOperator(int gridpoints, fftw_complex *op, const std::function<double(double)>& potential,
                             double space_width, double time_width, double start) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = cos(potential(i * space_width + start) * time_width / 2.0);
        op[i][1] = -sin(potential(i * space_width + start) * time_width / 2.0);
    }
}

// potential functions
// step potential
[[maybe_unused]] inline std::function<double(double)> step(double pos, double strength_1, double strength_2) {
    return [pos, strength_1, strength_2](double x) {
        return (x <= pos) ? strength_1 : strength_2;
    };
}
// harmonic oscillator
[[maybe_unused]] inline std::function<double(double)> sho(double strength) {
    return [strength](double x) {
        return strength*x*x/2;
    };
}
// well
[[maybe_unused]] inline std::function<double(double)> barrier(double start_pos, double end_pos, double strength) {
    return [start_pos, end_pos, strength](double x) {
        return (x <= end_pos && x >= start_pos) ? strength : 0;
    };
}
// triangle
[[maybe_unused]] inline std::function<double(double)> triangle(double strength) {
    return [strength](double x) {
        return (x > -3 && x <= 0) ? -strength*x : (x > 0 && x <= 3 ) ? strength*x : 0;
    };
}
// wall
[[maybe_unused]] inline std::function<double(double)> wall(double pos, double strength) {
    return [pos, strength](double x) {
        return (x <= pos) ? 0 : strength;
    };
}

#endif //PROPAGATE_POTENTIALS_H