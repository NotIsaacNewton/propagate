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

// real time propagation
// creates potential operator array from data file and outputs to op
void definePotentialOperator(int gridpoints, fftw_complex *op, const std::string& potfile, double time_width);

// calculates free-particle operator based on general values and outputs to op
void defineKineticOperator(int gridpoints, fftw_complex *op, double space_width, double time_width);

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

// imaginary time propagation
// calculates free-particle operator based on general values and outputs to op
void defineImKineticOperator(int gridpoints, fftw_complex *op, double space_width, double time_width);

// creates potential operator array from potential function and outputs to op
void defineImPotentialOperator(int gridpoints, fftw_complex *op,
    const std::function<double(double)>& potential, double space_width,
    double time_width, double start);

// propagates wavefunction based on general values
void ipropagate(const inputs& in, fftw_complex *psi, const std::function<double(double)>& potential, std::string& data);

#endif //PROPAGATE_H
