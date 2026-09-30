//
// Created by Arian Dovald on 9/18/25.
//

#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <cmath>
#include "fftw_complex_tools.h"

// scales entire array by a scalar
void scale_fftw_complex(const double scalar, fftw_complex *complex_vec, const int size) {
    for (int i = 0; i < size; i++) {
        complex_vec[i][0] *= scalar;
        complex_vec[i][1] *= scalar;
    }
}

// writes fftw_complex array to file
void fftw_complex_array_to_file(const double& start, const int& size, const double& width,
    const std::string& file, const fftw_complex *function) {
    std::ofstream potwrite;
    potwrite.open(file);
    if (potwrite.is_open()) {
        for (int i=0; i<size; i++) {
            potwrite << i*width + start << " " << function[i][0] << " " << function[i][1] << "\n";
        }
        potwrite.close();
    } else {
        std::cerr << "Failed to open " << file << ".\n";
    }
}

void fftw_complex_func_to_array(const double& start, const int& size, const double& width,
    const std::function<void(double, fftw_complex)>& function, fftw_complex *out) {
    fftw_complex temp;
    for (int i = 0; i < size; i++) {
        function(i*width + start, temp);
        out[i][0] = temp[0];
        out[i][1] = temp[1];
    }
}

// writes fftw_complex function to file
void fftw_complex_func_to_file(const inputs& in, const std::string& savefile,
    const std::function<void(double, fftw_complex)>& wavefunction) {
    // allocate temp array with RAII
    const auto temp = fftw_alloc_complex(in.space_grid);
    std::unique_ptr<fftw_complex, void(*)(void*)> psip{temp, fftw_free};
    // write wavefunction to array, save array to file
    fftw_complex_func_to_array(in.initial_pos,in.space_grid,in.dx,
        wavefunction, temp);
    fftw_complex_array_to_file(in.initial_pos, in.space_grid, in.dx,
        savefile, temp);
}

// writes array of fftw_complex functions to file
void fftw_complex_func_array_to_file(const inputs& in, const std::string& save_dir, const std::string& savefile,
    const std::vector<std::function<void(double, fftw_complex)>>& wavefunction) {
    for (int c = 0; c < in.channels; c++) {
        fftw_complex_func_to_file(in, save_dir + "/" + (savefile + "_") + std::to_string(c) + ".dat",
            wavefunction[c]);
    }
}

// reads to fftw_complex array from file
void fftw_complex_array_from_file(const std::string& file, fftw_complex *function, const int& size) {
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        std::string line;
        int n = 0;
        while (std::getline(read, line) && n < size) {
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
void print_fftw_complex(const int size, const fftw_complex *in) {
    for (int i = 0; i < size; i++) {
        std::cout << "(" << in[i][0] << ", " << in[i][1] << ")\n";
    }
}

// finds square amplitude of fftw_complex
void fftw_complex_square(const fftw_complex* function, std::vector<double>& out) {
    const std::size_t size = out.size();
    for (int i = 0; i < size; i++) {
        out[i] = function[i][0]*function[i][0] + function[i][1]*function[i][1];
    }
}

// integrates through array
double fftw_complex_integrate(const int size, const double width, const std::vector<double>& in) {
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += in[i]*width;
    }
    return sum;
}

// gets the norm of an fftw_complex vector
double norm(const int gridpoints, const double gridwidth, const fftw_complex* psi) {
    std::vector<double> psi_squared(gridpoints); // stores |psi|^2
    fftw_complex_square(psi, psi_squared); // calculates |psi|^2
    const double norm = fftw_complex_integrate(gridpoints, gridwidth, psi_squared); // calculate and store norm
    return norm;
}

// normalizes an fftw_complex vector
void normalize(const int gridpoints, const double gridwidth, fftw_complex* psi) {
    const double mag = norm(gridpoints, gridwidth, psi); // get norm
    scale_fftw_complex(1/sqrt(mag), psi, gridpoints); // normalize psi
}

// gets the norm of a vector of fftw_complex vectors, in parallel over the spatial grid
double normCC(const int gridpoints, const double dx, const int channels, const std::vector<fftw_complex*>& psi) {
    double total = 0.0;
    #pragma omp parallel for reduction(+:total)
    for (int i = 0; i < gridpoints; i++) {
        for (int c = 0; c < channels; c++) {
            total += psi[c][i][0]*psi[c][i][0] + psi[c][i][1]*psi[c][i][1];
        }
    }
    return total * dx;
}

// normalizes a vector of fftw_complex vectors, in parallel over channels
void normalizeCC(const int gridpoints, const double dx, const int channels, const std::vector<fftw_complex*>& psi) {
    const double mag = normCC(gridpoints, dx, channels, psi);
    const double scale = 1.0 / sqrt(mag);
    #pragma omp parallel for
    for (int c = 0; c < channels; c++) {
        scale_fftw_complex(scale, psi[c], gridpoints);
    }
}

// fftw_complex to std::complex<double>
std::complex<double> fftw_complex_to_std_complex(const fftw_complex& fftw) {
    return {fftw[0], fftw[1]};
}

// std::complex<double> to fftw_complex
void std_complex_to_fftw_complex(const std::complex<double>& std, fftw_complex& fftw) {
    fftw[0] = std.real();
    fftw[1] = std.imag();
}

// multiplies fftw_complex and std::complex<double>
std::complex<double> operator*(const double(&fftw)[2], const std::complex<double>& std) {
    return fftw_complex_to_std_complex(fftw) * std;
}
std::complex<double> operator*(const std::complex<double>& std, const double(&fftw)[2]) {
    return fftw_complex_to_std_complex(fftw) * std;
}

// adds two fftw_complex
std::complex<double> fftw_complex_add(const fftw_complex& a, const fftw_complex& b) {
    return {a[0] + b[0], a[1] + b[1]};
}
