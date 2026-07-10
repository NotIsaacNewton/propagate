//
// Created by Arian Dovald on 6/23/25.
//

#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <functional>
#include "filetools.h"
#include "fftw3.h"
#include "fftw_complex_tools.h"

// calculates momentum-squared grid based on general values
std::vector<double> psquared(int gridpoints, double space_width);

// creates potential operator array from data file and outputs to op
void definePotentialOperator(const inputs& in, fftw_complex *op, const std::string& potfile, bool imProp);

// calculates free-particle operator based on general values and outputs to op
void defineKineticOperator(const inputs& in, fftw_complex *op, bool imProp);

// potential energy exponential operator
void applyPotentialOperator(int gridpoints, fftw_complex *psi, const fftw_complex *V);

// kinetic energy exponential operator
void applyKineticOperator(int gridpoints, fftw_complex *psi, const fftw_complex *T);

// output to file
void writeOutput(const fftw_complex *psi, int t, int gridpoints, int output_ndx, int output_ndt,
    std::vector<double> &buffer);

// struct bundling all resources needed for using fftw in propagation loops
struct fftwResources {
    // pointers used as wrappers for RAII
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> fft_ptr;
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)> ifft_ptr;
    std::unique_ptr<fftw_complex, void(*)(void*)> Vp;
    std::unique_ptr<fftw_complex, void(*)(void*)> Tp;
};

fftwResources fftwPrep(const inputs& in, fftw_complex *psi, const std::string& data, bool imProp);

void propTick(int gridpoints, fftw_complex *psi, const fftw_complex* V, const fftw_complex* T,
    fftw_plan fft, fftw_plan ifft, double scale);

// propagates wavefunction based on general values
void propagate(const inputs& in, fftw_complex *psi, const std::string& data, bool imProp);

#endif //PROPAGATE_H
