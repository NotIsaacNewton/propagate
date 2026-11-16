//
// Created by Arian Dovald on 11/16/26.
//

#ifndef PROPAGATE_PROPAGATE_TD_CC_SEQ_H
#define PROPAGATE_PROPAGATE_TD_CC_SEQ_H

#include "propagate_td_cc.h"

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

// applies potential operator in parallel threads across spatial gridpoints
void applyPotentialOperatorCC(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi);
void applyPotentialOperatorCCReverse(const inputs& in, const std::vector<fftw_complex*>& diag,
    const std::vector<std::vector<std::vector<coupling>>>& coup, const std::vector<fftw_complex*>& psi);

// propagates psi in potential from tick to tick + 1
void propTickTDCCSeq(const int& tick, const inputs& in, const std::vector<fftw_complex*>& psi,
    const potCC& potential, const std::vector<fftw_complex*>& V_d,
    std::vector<std::vector<std::vector<coupling>>>& V_c, const fftw_complex* T,
    const std::vector<fftw_plan>& fft, const std::vector<fftw_plan>& ifft, double scale);

// sets up and propagates psi (here many channels) in a TD CC potential based on general values
void propagateTDCCSeq(const inputs& in, const std::string& data);

#endif //PROPAGATE_PROPAGATE_TD_CC_SEQ_H
