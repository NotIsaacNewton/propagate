//
// Created by Arian Dovald on 6/26/25.
//

#ifndef FFTW_COMPLEX_TOOLS_H
#define FFTW_COMPLEX_TOOLS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "armpl.h"
#include "fftw3.h"

// writes fftw_complex to file
void fftw_complex_array_to_file(const double& start, const double& end, const double& width, const std::string& file, fftw_complex *function) {
    const double n = (end-start)/width;
    std::ofstream potwrite;
    potwrite.open(file);
    if (potwrite.is_open()) {
        std::print("Writing to {}...\n",file);
        for (int i=0; i<n; i++) {
            potwrite << i*width + start << " " << function[i][0] << " " << function[i][1] << std::endl;
        }
        potwrite.close();
    } else {
        std::cerr << "Failed to open " << file << "." << std::endl;
    }
}

void fftw_complex_func_to_array(const double& start, const double& end, const double& width, void function(double x, fftw_complex temp), fftw_complex *out) {
    fftw_complex temp;
    const double n = (end-start)/width;
    for (int i = 0; i < n; i++) {
        function(i*width + start, temp);
        out[i][0] = temp[0];
        out[i][1] = temp[1];
    }
}

// reads to fftw_complex array from file
[[maybe_unused]] void fftw_complex_array_from_file(const std::string& file, fftw_complex *function) {
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        std::print("Reading from {}...\n",file);
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
        std::cerr << "Failed to open " << file << "." << std::endl;
    }
}

// prints fftw_complex (mostly for debugging)
[[maybe_unused]] void print_fftw_complex(int size, fftw_complex *in) {
    for (int i = 0; i < size; i++) {
        std::cout << "(" << in[i][0] << ", " << in[i][1] << ")" << std::endl;
    }
}

// finds square amplitude of fftw_complex
void fftw_complex_square(const int size, fftw_complex *function, double *out) {
    for (int i = 0; i < size; i++) {
        out[i] = function[i][0]*function[i][0] + function[i][1]*function[i][1];
    }
}

// integrates through array
double fftw_complex_integrate(const int size, double width, const double *in) {
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += in[i]*width;
    }
    return sum;
}

#endif //FFTW_COMPLEX_TOOLS_H
