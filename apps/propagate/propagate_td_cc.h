//
// Created by Arian Dovald on 9/2/26.
//

#ifndef PROPAGATE_PROPAGATE_TD_CC_H
#define PROPAGATE_PROPAGATE_TD_CC_H

#include "propagate.h"
#include "matrix_tools.h"

// gets potential and returns array of hermitian matrices
// first index is time*channel and second index is position*channel
std::vector<std::vector<hermitian_matrix>> getPotentialCC(const inputs& in, const std::string& data);

// creates array of diagonal potential operator arrays from potential at tick and outputs to op
// op at first level is indexed for each channel and at second level for spatial grid
void definePotentialOperatorDiag(const inputs& in, const int& tick, const std::vector<fftw_complex*>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential);
// creates array of diagonal potential operator arrays from potential at tick and outputs to op
void definePotentialOperatorDiag(const inputs& in, const int& tick, const std::vector<fftw_complex*>& op,
    const std::function<hermitian_matrix(double, double)>& potential);

// struct used to carry usable information of coupling elements (real-space two method)
// contains cos(|V_{ij}|*dt/2), sin(|V_{ij}|*dt/2), -i*V_{ij}/|V_{ij}|
struct coupling {
    std::complex<double> cos_factor;
    std::complex<double> sin_factor;
    std::complex<double> exp_phase;
};

// executes fft plans in parallel for all channels
void fftExecuteCC(int gridpoints, int channels, const std::vector<fftw_complex*>& psi,
    const std::vector<fftw_plan>& fft_plans);
// executes ifft plans in parallel for all channels
void ifftExecuteCC(int gridpoints, int channels, const std::vector<fftw_complex*>& psi,
    const std::vector<fftw_plan>& ifft_plans);

// applies diagonal part of potential operator
void applyDiagonalOperator(int point, int channels, const std::vector<fftw_complex*>& diag,
    const std::vector<fftw_complex*>& psi);

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

// variant for handling potential input types
using potCC = std::variant<std::vector<std::vector<hermitian_matrix>>,
    std::function<hermitian_matrix(double, double)>>;

// gets initial wavepacket handled with RAII, outputing a vector of unique pointers
std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> getWavepacketCC(const inputs& in, const std::string& data);

// writes output serially across channels (cannot be done in parallel)
void writeOutputCC(const std::vector<fftw_complex*>& psi, int t, const inputs& in, std::vector<double>& buffer);

// calculate and save PEY to file
void savePEY(const std::vector<fftw_complex*>& psi, const inputs& in);

#endif //PROPAGATE_PROPAGATE_TD_CC_H
