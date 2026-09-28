//
// Created by Arian Dovald on 6/26/25.
//

#ifndef FFTW_COMPLEX_TOOLS_H
#define FFTW_COMPLEX_TOOLS_H

#include <string>
#include "fftw3.h"
#include "file_tools.h"
#include <complex>

// scales entire array by a scalar
void scale_fftw_complex(double scalar, fftw_complex *complex_vec, int size);

// writes fftw_complex array to file
void fftw_complex_array_to_file(const double& start, const int& size, const double& width,
    const std::string& file, const fftw_complex *function);

void fftw_complex_func_to_array(const double& start, const int& size, const double& width,
    const std::function<void(double, fftw_complex)>& function, fftw_complex *out);

// writes fftw_complex function to file
void fftw_complex_func_to_file(const inputs& in, const std::string& savefile,
    const std::function<void(double, fftw_complex)>& wavefunction);

// reads to fftw_complex array from file
void fftw_complex_array_from_file(const std::string& file, fftw_complex *function, const int& size);

// prints fftw_complex (mostly for debugging)
void print_fftw_complex(int size, const fftw_complex *in);

// finds square amplitude of fftw_complex
void fftw_complex_square(const fftw_complex* function, std::vector<double>& out);

// integrates through array
double fftw_complex_integrate(int size, double width, const std::vector<double>& in);

// gets the norm of an fftw_complex vector
double norm(int gridpoints, double gridwidth, const fftw_complex* psi);

// normalizes an fftw_complex vector
void normalize(int gridpoints, double gridwidth, fftw_complex* psi);

// gets the norm of a vector of fftw_complex vectors, in parallel over the spatial grid
double normCC(int gridpoints, double dx, int channels, const std::vector<fftw_complex*>& psi);

// normalizes a vector of fftw_complex vectors, in parallel over channels
void normalizeCC(int gridpoints, double dx, int channels, const std::vector<fftw_complex*>& psi);

// fftw_complex to std::complex<double>
std::complex<double> fftw_complex_to_std_complex(const fftw_complex& fftw);

// std::complex<double> to fftw_complex
void std_complex_to_fftw_complex(const std::complex<double>& std, fftw_complex& fftw);

// multiplies fftw_complex and std::complex<double>
std::complex<double> operator*(const double(&fftw)[2], const std::complex<double>& std);
std::complex<double> operator*(const std::complex<double>& std, const double(&fftw)[2]);

// adds two fftw_complex
std::complex<double> fftw_complex_add(const fftw_complex& a, const fftw_complex& b);

#endif //FFTW_COMPLEX_TOOLS_H
