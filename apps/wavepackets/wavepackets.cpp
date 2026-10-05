//
// Created by Arian Dovald on 9/9/25.
//

#include <functional>
#include <numbers>
#include <string>
#include <cmath>
#include <map>
#include <unordered_map>
#include "wavepackets.h"
#include "file_tools.h"
#include "interpolate_1d.h"

double gaussian(const double x, const double delta, const double pos) {
    return pow(1/(std::numbers::pi * (delta*delta)), 1.0/4.0)*exp(-(x-pos)*(x-pos)/(2*(delta*delta)));
}

// gives Gaussian momentum
void gaussianMoving(const double x, const double delta, const double momentum, const double pos,
    fftw_complex out) {
    const double re = cos(momentum * x);
    const double im = sin(momentum * x);
    out[0] = re*gaussian(x, delta, pos);
    out[1] = im*gaussian(x, delta, pos);
}

// prepares wave-packet
std::function<void(double, fftw_complex)> gaussianWP(const double delta, const double momentum,
    const double pos) {
    return [delta, momentum, pos](const double x, fftw_complex out) {
        gaussianMoving(x, delta, momentum, pos, out);
    };
}

// SHO ground-state
std::function<void(double, fftw_complex)> shoGround() {
    return [](const double x, fftw_complex out) {
        const double psi = pow(1/std::numbers::pi, 1.0/4.0)*exp(-x*x/2);
        out[0] = psi;
        out[1] = 0;
    };
}

// SHO first excited state
std::function<void(double, fftw_complex)> shoExcited() {
    return [](const double x, fftw_complex out) {
        const double psi = pow(1/std::numbers::pi, 1.0/4.0)*sqrt(2)*x*exp(-x*x/2);
        out[0] = psi;
        out[1] = 0;
    };
}

// test state (for imaginary time propagation)
std::function<void(double, fftw_complex)> test() {
    return [](const double x, fftw_complex out) {
        double psi;
        x <= 5 && x >= -5 ? psi = 1 : psi = 0;
        out[0] = psi;
        out[1] = 0;
    };
}

// utility function for filling intially-empty channels
std::function<void(double, fftw_complex)> zeroState() {
    return [](const double, fftw_complex out) {
        out[0] = 0.0;
        out[1] = 0.0;
    };
}

// read function from file of raw doubles
std::function<void(double, fftw_complex)> psiFromFile(const inputs& in) {
    std::vector<double> temp(in.space_grid_coarse); // stores temp wavefunction for interpolation
    readArray1D(in.input_psi_file, temp); // reads from file
    std::vector<double> grid(in.space_grid_coarse); // stores coarse grid on which wavefunction is defined
    const double dx = (in.final_pos-in.initial_pos)/(in.space_grid_coarse-1); // width of coarse grid
    // write grid
    #pragma omp parallel for default(none) shared(in, grid, dx)
    for (int i = 0; i < in.space_grid_coarse; i++) {
        grid[i] = in.initial_pos + i*dx;
    }
    spline_interp interpolator(grid, temp); // spline interpolation object
    return [interpolator](const double x, fftw_complex out) {
        out[0] = interpolator.interp(x);
        out[1] = 0.0;
    };
}

// map of options
std::function<void(double, fftw_complex)> buildWavepacket(const inputs& in) {
    const std::unordered_map<std::string, std::function<std::function<void(double, fftw_complex)>()>> wavepackets = {
        {"gaussian",    [in] {
            return gaussianWP(in.delta,in.momentum, in.wp_position);
            }},
        {"shoground",   [] { return shoGround(); }},
        {"shoexcited",  [] { return shoExcited(); }},
        {"test",        [] { return test(); }},
        {"file",        [in] { return psiFromFile(in); }},
    };
    return wavepackets.at(in.wavepacket_type)();
}
