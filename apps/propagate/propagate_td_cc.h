//
// Created by Arian Dovald on 9/2/26.
//

#ifndef PROPAGATE_PROPAGATE_TD_CC_H
#define PROPAGATE_PROPAGATE_TD_CC_H

#include "propagate.h"
#include "matrix_tools.h"

// gets potential and returns array of hermitian matrices
// NOTE: first index is time*channel and second index is position*channel
std::vector<std::vector<hermitian_matrix>> getPotentialCC(const inputs& in, const std::string& data);

// creates array of diagonal potential operator arrays from potential at tick and outputs to op
// NOTE: op at first level is indexed for each channel and at second level for spatial grid
void definePotentialOperatorDiag(const inputs& in, const int& tick, const std::vector<fftw_complex*>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential);

// struct used to carry usable information of coupling elements
// NOTE: should contain cos(|V_{ij}|*dt/2), sin(|V_{ij}|*dt/2), -i*V_{ij}/|V_{ij}|
struct coupling {
    std::complex<double> cos_factor;
    std::complex<double> sin_factor;
    std::complex<double> exp_phase;
};

// creates array of coupling operator arrays from data at tick and outputs to op
// NOTE: op at first level is indexed for spatial grid and at second and third levels for each pair of channels
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<std::vector<coupling>>>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential);

// applies coupling part of potential operator for all channels at a specific gridpoint
void applyCouplingOperator(int point, int channels, const std::vector<std::vector<std::vector<coupling>>>& coup,
    const std::vector<fftw_complex*>& psi);

void applyDiagonalOperator(int point, int channels, const std::vector<fftw_complex*>& diag,
    const std::vector<fftw_complex*>& psi);

// applies potential operator in parallel threads across spatial gridpoints
void applyPotentialOperatorCC(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi);
void applyPotentialOperatorCCReverse(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi);

// executes fft plans in parallel for all channels
void fftExecuteCC(int gridpoints, int channels, const std::vector<fftw_complex*>& psi,
    const std::vector<fftw_plan>& fft_plans);

// executes ifft plans in parallel for all channels
void ifftExecuteCC(int gridpoints, int channels, const std::vector<fftw_complex*>& psi,
    const std::vector<fftw_plan>& ifft_plans);

// applies kinetic energy operator using applyKineticOperator from propagate.cpp in parallel
void applyKineticOperatorCC(int gridpoints, int channels, const std::vector<fftw_complex*>& psi, const fftw_complex *T);

// bundles all RAII-sensitive resources for propagation and fftw
struct fftwResourcesCC {
    std::vector<std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>> fft_ptrs;
    std::vector<std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>> ifft_ptrs;
    std::unique_ptr<fftw_complex, void(*)(void*)> Tp;
    std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> V_d;
    std::vector<std::vector<std::vector<coupling>>> V_c;
};

// prepares fftw and propagation variables with RAII
fftwResourcesCC fftwPrepTDCC(const inputs& in, const std::vector<fftw_complex*>& psi, const std::string& data);

// normalizes fftw results in parallel for all channels
void fftwNormCC(const inputs& in, double scale, const std::vector<fftw_complex*>& psi);

// propagates psi in potential from tick to tick + 1
void propTickTDCC(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const std::vector<std::vector<hermitian_matrix>>& potential, const std::vector<fftw_complex*>& V_d,
    std::vector<std::vector<std::vector<coupling>>>& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale);

// gets initial wavepacket handled with RAII, outputing a vector of unique pointers
std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> getWavepacketCC(const inputs& in, const std::string& data);

// writes output serially across channels (cannot be done in parallel)
void writeOutputCC(const std::vector<fftw_complex*>& psi, int t, const inputs& in, std::vector<double>& buffer);

// sets up and propagates psi (here many channels) in a TD CC potential based on general values
void propagateTDCC(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPAGATE_TD_CC_H
