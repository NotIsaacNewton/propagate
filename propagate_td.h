//
// Created by Arian Dovald on 8/19/26.
//

#ifndef PROPAGATE_PROPOGATE_TD_H
#define PROPAGATE_PROPOGATE_TD_H

#include "interpolate_1d.h"
#include "propagate.h"

// creates potential operator array from data file and outputs to op
void definePotentialOperatorTD(const inputs& in, fftw_complex *op, const std::vector<double>& potential);

// prepares fftw variables with RAII
fftwResources fftwPrepTD(const inputs& in, fftw_complex *psi, const std::string& data);

struct interpolator {
    std::vector<double> pos_grid;
    std::vector<spline_interp> time_interp;
};

// prepares interpolator for potential
interpolator prepInterpolation(const inputs& in, const std::vector<std::vector<double>>& potential);

void interpolatePotential(const inputs& in, interpolator& interp, std::vector<double>& pot_coarse,
    std::vector<double>& potential, int time);

void propTickTD(int tick, const inputs& in, fftw_complex *psi, std::vector<double>& pot_coarse, std::vector<double>& potential,
    fftw_complex* V, const fftw_complex* T, fftw_plan fft, fftw_plan ifft, interpolator& interp, double scale);

void propagateTD(const inputs& in, fftw_complex *psi, const std::string& data);

#endif //PROPAGATE_PROPOGATE_TD_H
