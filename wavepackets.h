//
// Created by Arian Dovald on 6/30/25.
//

#ifndef WAVEPACKETS_H
#define WAVEPACKETS_H

#include <functional>
#include <numbers>
#include "fftw3.h"

// Gaussian wave-packet
// Gaussian
[[maybe_unused]] inline double gaussian(const double x, const double delta, const double pos) {
    return pow(1/(std::numbers::pi * (delta*delta)), 1.0/4.0)*exp(-(x-pos)*(x-pos)/(2*(delta*delta)));
}
// gives Gaussian momentum
[[maybe_unused]] inline void gaussianMoving(const double x, const double delta, const double momentum, const double pos,
    fftw_complex out) {
    const double re = cos(momentum * x);
    const double im = sin(momentum * x);
    out[0] = re*gaussian(x, delta, pos);
    out[1] = im*gaussian(x, delta, pos);
}
// prepares wave-packet
[[maybe_unused]] inline std::function<void(double, fftw_complex)> gaussianWP(const double delta, const double momentum,
    const double pos) {
    return [delta, momentum, pos](const double x, fftw_complex out) {
        gaussianMoving(x, delta, momentum, pos, out);
    };
}

// SHO ground-state
[[maybe_unused]] inline std::function<void(double, fftw_complex)> shoGround() {
    return [](const double x, fftw_complex out) {
        const double psi = pow(1/std::numbers::pi, 1.0/4.0)*exp(-x*x/2);
        out[0] = psi;
        out[1] = 0;
    };
}

// SHO first excited state
[[maybe_unused]] inline std::function<void(double, fftw_complex)> shoExcited() {
    return [](const double x, fftw_complex out) {
        const double psi = pow(1/std::numbers::pi, 1.0/4.0)*sqrt(2)*x*exp(-x*x/2);
        out[0] = psi;
        out[1] = 0;
    };
}

// test state (for imaginary time propagation)
[[maybe_unused]] inline std::function<void(double, fftw_complex)> test() {
    return [](const double x, fftw_complex out) {
        double psi;
        (x <= 5 && x >= -5) ? psi = 1 : psi = 0;
        out[0] = psi;
        out[1] = 0;
    };
}

#endif //WAVEPACKETS_H
