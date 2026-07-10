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
void writeOutput(const fftw_complex *psi, int t, double start, int gridpoints,
    double space_width, double time_width, int output_ndx, int output_ndt,
    std::ostringstream &buffer);

// propagates wavefunction based on general values
void propagate(const inputs& in, fftw_complex *psi, const std::string& data);

#endif //PROPAGATE_H
