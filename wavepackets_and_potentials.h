//
// Created by Arian Dovald on 6/30/25.
//

#ifndef WAVEPACKETS_AND_POTENTIALS_H
#define WAVEPACKETS_AND_POTENTIALS_H

#include <functional>
#include "armpl.h"
#include "fftw3.h"

// Gaussian wave-packet
// Gaussian
[[maybe_unused]] double gaussian(double x, const double delta) {
    return pow(1/(M_PI * (delta*delta)), 1.0/4.0)*exp(-(x+3)*(x+3)/(2*(delta*delta)));
}
// gives Gaussian momentum
[[maybe_unused]] void gaussianMoving(double x, const double delta, const double momentum, fftw_complex out) {
    double re = cos(momentum * x);
    double im = sin(momentum * x);
    out[0] = re*gaussian(x, delta);
    out[1] = im*gaussian(x, delta);
}
// prepares wave-packet (narrows parameters because the function from fftw_complex_tools.h is too picky)
// TODO: make fftw_complex_tools.h less picky
[[maybe_unused]] void gaussianWP(double x, fftw_complex out) {
    gaussianMoving(x, 0.5, 15, out);
}

// SHO ground-state
[[maybe_unused]] void shoGround(double x, fftw_complex out) {
    double psi = pow(1/M_PI, 1.0/4.0)*exp(-x*x/2);
    out[0] = psi;
    out[1] = 0;
}
// SHO fist excited state
[[maybe_unused]] void shoExcited(double x, fftw_complex out) {
    double psi = pow(1/M_PI, 1.0/4.0)*sqrt(2)*x*exp(-x*x/2);
    out[0] = psi;
    out[1] = 0;
}

// test state (for imaginary time propagation)
[[maybe_unused]] void test(double x, fftw_complex out) {
    double psi;
    (x <= 5 && x >= -5) ? psi = 1 : psi = 0;
    out[0] = psi;
    out[1] = 0;
}

// potential functions
// step potential
[[maybe_unused]] std::function<double(double)> step(double pos, double strength_1, double strength_2) {
    return [pos, strength_1, strength_2](double x) {
        return (x <= pos) ? strength_1 : strength_2;
    };
}
// harmonic oscillator
[[maybe_unused]] double sho(double x) {
    return x*x/2;
}
// well
[[maybe_unused]] std::function<double(double)> barrier(double start_pos, double end_pos, double strength) {
    return [start_pos, end_pos, strength](double x) {
        return (x <= end_pos && x >= start_pos) ? strength : 0;
    };
}
// triangle
[[maybe_unused]] double triangle(double x) {
    return (x > -3 && x <= 0) ? -x : (x > 0 && x <= 3 ) ? x : 0;
}
// wall
[[maybe_unused]] std::function<double(double)> wall(double pos, double strength) {
    return [pos, strength](double x) {
        return (x <= pos) ? 0 : strength;
    };
}

#endif //WAVEPACKETS_AND_POTENTIALS_H
