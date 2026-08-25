//
// Created by Arian Dovald on 8/19/26.
//

#ifndef PROPAGATE_PROPOGATE_TD_H
#define PROPAGATE_PROPOGATE_TD_H

#include "interpolate_1d.h"
#include "propagate.h"

// gets potential and returns arrays
std::vector<std::vector<double>> getPotential(const inputs& in, const std::string& data);

// creates potential operator array from data file and outputs to op
void definePotentialOperatorTD(const inputs& in, const int& tick, fftw_complex *op,
    const std::vector<std::vector<double>>& potential);

// prepares fftw variables with RAII
fftwResources fftwPrepTD(const inputs& in, fftw_complex *psi, const std::string& data);

// propagates psi in potential from tick to tick + 1
void propTickTD(const int& tick, const inputs& in, fftw_complex *psi, const std::vector<std::vector<double>>& potential,
    fftw_complex* V, const fftw_complex* T, const fftw_plan* fft, const fftw_plan* ifft, double scale);

// sets up and propagates psi in a TD potential based on general values
void propagateTD(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPOGATE_TD_H
