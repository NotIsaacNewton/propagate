//
// Created by Arian Dovald on 6/26/25.
//

#ifndef FFTW_COMPLEX_TOOLS_H
#define FFTW_COMPLEX_TOOLS_H

#include <string>
#include "fftw3.h"
#include "filetools.h"

// scales entire array by a scalar
void scale_fftw_complex(double scalar, fftw_complex *complex_vec, int size);

// writes fftw_complex array to file
void fftw_complex_array_to_file(const double& start, const double& end, const double& width,
    const std::string& file, const fftw_complex *function);

void fftw_complex_func_to_array(const double& start, const double& end, const double& width,
    const std::function<void(double, fftw_complex)>& function, fftw_complex *out);

// writes fftw_complex function to file
void fftw_complex_func_to_file(const inputs& in, const std::string& savefile,
    const std::function<void(double, fftw_complex)>& wavefunction);

// reads to fftw_complex array from file
void fftw_complex_array_from_file(const std::string& file, fftw_complex *function);

// prints fftw_complex (mostly for debugging)
void print_fftw_complex(int size, const fftw_complex *in);

// finds square amplitude of fftw_complex
void fftw_complex_square(const fftw_complex* function, std::vector<double>& out);

// integrates through array
double fftw_complex_integrate(int size, double width, const std::vector<double>& in);

#endif //FFTW_COMPLEX_TOOLS_H
