//
// Created by Arian Dovald on 9/2/26.
//

#ifndef PROPAGATE_PROPAGATE_TD_CC_H
#define PROPAGATE_PROPAGATE_TD_CC_H

#include "propagate.h"
#include "matrix_tools.h"

// gets potential and returns array of hermitian matrices
// NOTE: first index is time and second index is position
//  files will follow naming convention potential_i_j.dat, allowing for parallel reads and writes, and reusability of
//  writeFunction2D in potential_td_cc.cpp and readArray2D in here
std::vector<std::vector<hermitian_matrix>> getPotentialCC(const inputs& in, const std::string& data);

// creates array potential operator arrays from data and outputs to op
// NOTE: op at first level is indexed for each channel and at second level for spatial grid
void definePotentialOperatorDiag(const inputs& in, const int& tick, std::vector<fftw_complex*>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential);

// struct used to carry usable information of coupling elements
// NOTE: should contain cos(|V_{ij}|*dt/2), sin(|V_{ij}|*dt/2), -i*V_{ij}/|V_{ij}|
struct coupling {
    double cos_factor;
    double sin_factor;
    std::complex<double> exp_phase;
};

// creates array of coupling operator arrays from data and outputs to op
// NOTE: op at first level is indexed for spatial grid and at second level is indexed for each pair of channels
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<coupling>>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential);

// applies potential operators across all channels in parallel, using applyPotentialOperator from propagate.cpp
void applyPotentialOperatorDiag(const inputs& in, std::vector<fftw_complex*>& potential,
    std::vector<fftw_complex*>& psi);

// applies coupling part of potential operator sequentially by pairs, and parallel across spatial gridpoints
// TODO: define multiplication between std::complex<double> and fftw_complex in fftw_complex_tools.cpp/.h
void applyCouplingOperator(const inputs& in, std::vector<std::vector<coupling>>& coupling,
    std::vector<fftw_complex*>& psi);

// applies kinetic energy operator using applyKineticOperator from propagate.cpp in parallel
void applyKineticOperatorCC(int gridpoints, std::vector<fftw_complex*>& psi, const fftw_complex *T);

struct fftwResourcesCC {};

// prepares fftw variables with RAII
fftwResourcesCC fftwPrepTDCC(const inputs& in, std::vector<fftw_complex*>& psi, const std::string& data);

// propagates psi in potential from tick to tick + 1
void propTickTDCC(const int& tick, const inputs& in, std::vector<fftw_complex*>& psi,
    const std::vector<std::vector<hermitian_matrix>>& potential,
    std::vector<fftw_complex*>& V_d, std::vector<std::vector<coupling>>& V_c, const fftw_complex* T,
    std::vector<const fftw_plan*> fft, std::vector<const fftw_plan*> ifft, double scale);

// sets up and propagates psi (here many channels) in a TD CC potential based on general values
void propagateTDCC(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPAGATE_TD_CC_H
