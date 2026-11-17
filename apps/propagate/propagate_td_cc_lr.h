//
// Created by Arian Dovald on 11/16/26.
//

#ifndef PROPAGATE_PROPAGATE_TD_CC_LR_H
#define PROPAGATE_PROPAGATE_TD_CC_LR_H

#include "propagate_td_cc.h"
#include <complex>

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

// propagates psi in potential from tick to tick + 1
void propTickLR(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const std::vector<fftw_complex*>& V_d,
    const coupling_lr& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale);

// sets up and propagates psi (here many channels) in a TD CC potential based on general values
void propagateTDCCLR(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPAGATE_TD_CC_LR_H
