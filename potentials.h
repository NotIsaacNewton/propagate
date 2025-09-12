//
// Created by Arian Dovald on 9/8/25.
//

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <functional>
#include "fftw3.h"

// creates potential operator array from potential function and outputs to op
[[maybe_unused]] inline void definePotentialOperator(const int gridpoints, fftw_complex *op,
    const std::function<double(double)>& potential, const double space_width,
    const double time_width, const double start) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = cos(potential(i * space_width + start) * time_width / 2.0);
        op[i][1] = -sin(potential(i * space_width + start) * time_width / 2.0);
    }
}

// potential functions
// step potential
[[maybe_unused]] inline std::function<double(double)> step(const double pos,
    const double strength_1, const double strength_2) {
    return [pos, strength_1, strength_2](const double x) {
        return (x <= pos) ? strength_1 : strength_2;
    };
}
// harmonic oscillator
[[maybe_unused]] inline std::function<double(double)> sho(const double strength) {
    return [strength](const double x) {
        return strength*x*x/2;
    };
}
// well
[[maybe_unused]] inline std::function<double(double)> barrier(const double start_pos,
    const double end_pos, const double strength) {
    return [start_pos, end_pos, strength](const double x) {
        return (x <= end_pos && x >= start_pos) ? strength : 0;
    };
}
// triangle
[[maybe_unused]] inline std::function<double(double)> triangle(const double strength) {
    return [strength](const double x) {
        return (x > -3 && x <= 0) ? -strength*x : (x > 0 && x <= 3 ) ? strength*x : 0;
    };
}
// wall
[[maybe_unused]] inline std::function<double(double)> wall(const double pos, const double strength) {
    return [pos, strength](const double x) {
        return (x <= pos) ? 0 : strength;
    };
}

#endif //POTENTIALS_H