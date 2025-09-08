//
// Created by Arian Dovald on 6/26/25.
//

#ifndef FFTW_COMPLEX_TOOLS_H
#define FFTW_COMPLEX_TOOLS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "fftw3.h"

// scales entire array by a scalar
inline void scale_fftw_complex(double scalar, fftw_complex *complex_vec, int size) {
    for (int i = 0; i < size; i++) {
        complex_vec[i][0] *= scalar;
        complex_vec[i][1] *= scalar;
    }
}

// writes fftw_complex to file
inline void fftw_complex_array_to_file(const double& start, const double& end, const double& width,
    const std::string& file, const fftw_complex *function) {
    const double n = (end-start)/width;
    std::ofstream potwrite;
    potwrite.open(file);
    if (potwrite.is_open()) {
        for (int i=0; i<n; i++) {
            potwrite << i*width + start << " " << function[i][0] << " " << function[i][1] << "\n";
        }
        potwrite.close();
    } else {
        std::cerr << "Failed to open " << file << ".\n";
    }
}

inline void fftw_complex_func_to_array(const double& start, const double& end, const double& width,
    const std::function<void(double, fftw_complex)>& function, fftw_complex *out) {
    fftw_complex temp;
    const double n = (end-start)/width;
    for (int i = 0; i < n; i++) {
        function(i*width + start, temp);
        out[i][0] = temp[0];
        out[i][1] = temp[1];
    }
}

// reads to fftw_complex array from file
[[maybe_unused]] inline void fftw_complex_array_from_file(const std::string& file, fftw_complex *function) {
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        std::string line;
        int n = 0;
        while (std::getline(read, line)) {
            std::istringstream readline(line);
            double x;
            double re;
            double im;
            readline >> x >> re >> im;
            function[n][0] = re;
            function[n][1] = im;
            ++n;
        }
        read.close();
    } else {
        std::cerr << "Failed to open " << file << ".\n";
        exit(1);
    }
}

// prints fftw_complex (mostly for debugging)
[[maybe_unused]] inline void print_fftw_complex(int size, const fftw_complex *in) {
    for (int i = 0; i < size; i++) {
        std::cout << "(" << in[i][0] << ", " << in[i][1] << ")\n";
    }
}

// finds square amplitude of fftw_complex
inline void fftw_complex_square(const int size, const fftw_complex *function, double *out) {
    for (int i = 0; i < size; i++) {
        out[i] = function[i][0]*function[i][0] + function[i][1]*function[i][1];
    }
}

// integrates through array
inline double fftw_complex_integrate(const int size, double width, const double *in) {
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += in[i]*width;
    }
    return sum;
}

#endif //FFTW_COMPLEX_TOOLS_H
