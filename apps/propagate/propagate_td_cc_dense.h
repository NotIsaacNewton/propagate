//
// Created by Arian Dovald on 11/16/26.
//

#ifndef PROPAGATE_PROPAGATE_TD_CC_DENSE_H
#define PROPAGATE_PROPAGATE_TD_CC_DENSE_H

#include "propagate_td_cc.h"

// applies potential operator using dense matrix diagonalization ZHEEVR
void applyPotentialOperatorDense(const inputs& in, int tick,
    const std::function<hermitian_matrix(double,double)>& potential, const std::vector<fftw_complex*>& psi);

// propagates psi in potential from tick to tick + 1
void propTickDense(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const std::function<hermitian_matrix(double,double)>& potential, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale);

// sets up and propagates psi (here many channels) in a TD CC potential based on general values
void propagateTDCCDense(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPAGATE_TD_CC_DENSE_H
