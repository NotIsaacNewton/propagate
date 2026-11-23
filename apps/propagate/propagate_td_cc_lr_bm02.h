//
// Created by Arian Dovald on 11/17/26.
//

#ifndef PROPAGATE_PROPAGATE_TD_CC_LR_BM02_H
#define PROPAGATE_PROPAGATE_TD_CC_LR_BM02_H

#include "propagate_td_cc_lr.h"

// apply low-rank coupling operator with modified coefficient alpha
void applyCouplingOperatorLRCoeff(const inputs& in, const coupling_lr& coup,
    const std::vector<fftw_complex*>& psi, int point, double tick, double alpha);

// applies potential operator for low-rank approximation with coefficient alpha and fractional tick
void applyPotentialOperatorLRCoeff(const inputs& in, double tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi, double alpha);
void applyPotentialOperatorLRCoeffReverse(const inputs& in, double tick, const std::vector<fftw_complex*>& diag,
    const coupling_lr& coup, const std::vector<fftw_complex*>& psi, double alpha);

// propagates psi in potential from tick to tick + 1
void propTickLRBM02(const double& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const std::function<hermitian_matrix(double,double)>& potential, const std::vector<fftw_complex*>& V_d,
    const coupling_lr& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale,
    const std::vector<double>& c);

// propagates wavefunction based on general values
void propagateTDCCLRBM02(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPAGATE_TD_CC_LR_BM02_H
