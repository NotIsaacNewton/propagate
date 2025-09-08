//
// Created by Arian Dovald on 9/8/25.
//

#ifndef PROPAGATE_POTENTIALS_H
#define PROPAGATE_POTENTIALS_H

#include <functional>

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