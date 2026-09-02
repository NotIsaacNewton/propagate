//
// Created by Arian Dovald on 6/23/25.
//

#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <functional>
#include <memory>
#include "../../lib/file_tools.h"
#include "fftw3.h"

// calculates momentum-squared grid based on general values
std::vector<double> psquared(int gridpoints, double space_width);

// gets wavepacket from file and puts it into a managed pointer
std::unique_ptr<fftw_complex, void(*)(void*)> getWavepacket(const inputs& in, const std::string& data);

// creates potential operator array from data file and outputs to op
void definePotentialOperator(const inputs& in, fftw_complex *op, const std::string& potfile, bool imProp);

// reshapes potential using interpolation
void reshapePotential(const inputs& in, std::vector<double> potential);

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

// prepares fftw resources needed for the propagation
fftwResources fftwPrep(const inputs& in, fftw_complex *psi, const std::string& data, bool imProp);

// propagates wavefunction for one tick in V
void propTick(int gridpoints, fftw_complex *psi, const fftw_complex* V, const fftw_complex* T,
    const fftw_plan* fft, const fftw_plan* ifft, double scale);

// propagates wavefunction based on general values
void propagate(const inputs& in, const std::string& data, bool imProp);

#endif //PROPAGATE_H
