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
// creates array of coupling operator arrays from data at tick and outputs to op
// op at first level is indexed for spatial grid and at second and third levels for each pair of channels
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<std::vector<coupling>>>& op,
    const std::vector<std::vector<hermitian_matrix>>& potential);
// creates array of coupling operator arrays from data at tick and outputs to op
void defineCouplingOperator(const inputs& in, const int& tick, std::vector<std::vector<std::vector<coupling>>>& op,
    const std::function<hermitian_matrix(double,double)>& potential);
// applies coupling part of potential operator for all channels at a specific gridpoint
void applyCouplingOperator(int point, int channels, const std::vector<std::vector<std::vector<coupling>>>& coup,
    const std::vector<fftw_complex*>& psi);
void applyCouplingOperatorReverse(int point, int channels, const std::vector<std::vector<std::vector<coupling>>>& coup,
    const std::vector<fftw_complex*>& psi);

// applies diagonal part of potential operator
void applyDiagonalOperator(int point, int channels, const std::vector<fftw_complex*>& diag,
    const std::vector<fftw_complex*>& psi);

// applies potential operator in parallel threads across spatial gridpoints
void applyPotentialOperatorCC(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi);
void applyPotentialOperatorCCReverse(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi);

// struct used to carry usable information of coupling elements (low-rank approximation)
struct coupling_lr {
    diag_matrix<double> diag_coup;
    std::vector<lapack_complex_double> eigenvectors;
};
// prep low-rank approximation method
coupling_lr prepCouplingOperatorLR(const inputs& in,
    const std::function<hermitian_matrix(double,double)>& potential);
coupling_lr prepCouplingOperatorLR(const inputs& in,
    const std::vector<std::vector<hermitian_matrix>>& potential);
// apply potential operator with low-rank approximation
void applyCouplingOperatorLR(const inputs& in, const coupling_lr& coup,
    const std::vector<fftw_complex*>& psi, int point, int tick);
// applies potential operator in parallel threads across spatial gridpoints for low-rank approximation
void applyPotentialOperatorLR(const inputs& in, int tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi);
// applies potential operator in parallel threads across spatial gridpoints for low-rank approximation
void applyPotentialOperatorLRReverse(const inputs& in, int tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi);

// applies potential operator using dense matrix diagonalization ZHEEVR
void applyPotentialOperatorDense(const inputs& in, int tick,
    const std::vector<std::vector<hermitian_matrix>>& potential, const std::vector<fftw_complex*>& psi);
void applyPotentialOperatorDense(const inputs& in, int tick,
    const std::function<hermitian_matrix(double,double)>& potential, const std::vector<fftw_complex*>& psi);

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
using potCC = std::variant<std::vector<std::vector<hermitian_matrix>>,
    std::function<hermitian_matrix(double, double)>>;
void propTickTDCC(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const std::vector<fftw_complex*>& V_d,
    std::vector<std::vector<std::vector<coupling>>>& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale);
void propTickDense(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale);
void propTickLR(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const std::vector<fftw_complex*>& V_d,
    const coupling_lr& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale);

// gets initial wavepacket handled with RAII, outputing a vector of unique pointers
std::vector<std::unique_ptr<fftw_complex, void(*)(void*)>> getWavepacketCC(const inputs& in, const std::string& data);

// writes output serially across channels (cannot be done in parallel)
void writeOutputCC(const std::vector<fftw_complex*>& psi, int t, const inputs& in, std::vector<double>& buffer);

// calculate and save PEY to file
void savePEY(const std::vector<fftw_complex*>& psi, const inputs& in);

// sets up and propagates psi (here many channels) in a TD CC potential based on general values
void propagateTDCC(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPAGATE_TD_CC_H
