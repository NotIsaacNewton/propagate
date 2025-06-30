//
// Created by Arian Dovald on 6/23/25.
//

#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include "armpl.h"
#include "fftw3.h"
#include "fftw_complex_tools.h"

// ANSI escape code for green
#define GREEN   "\033[32m"

// scales entire array by a scalar
void scale_fftw_complex(double scalar, fftw_complex *complex_vec, int size) {
    for (int i = 0; i < size; i++) {
        complex_vec[i][0] *= scalar;
        complex_vec[i][1] *= scalar;
    }
}

// calculates momentum-squared based on general values
double psquared(int i, int gridpoints, double space_width) {
    double p = 2 * M_PI * i / (gridpoints * space_width) - (M_PI / space_width);
    return p*p;
}

// real time propagation
// calculates free-particle operator based on general values and outputs to op
void defH(int gridpoints, fftw_complex *op, double space_width, double time_width) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = cos(time_width * psquared(i, gridpoints, space_width) / 2);
        op[i][1] = -sin(time_width * psquared(i, gridpoints, space_width) / 2);
    }
}

// creates potential operator array from potential function and outputs to op
void defV(int gridpoints, fftw_complex *op, double potential(double x), double space_width, double time_width, double start) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = cos(potential(i * space_width + start) * time_width / 2.0);
        op[i][1] = -sin(potential(i * space_width + start) * time_width / 2.0);
    }
}

// propogates wavefunction based on general values
// TODO: allow for time-dependent potentials (put in grid)
[[maybe_unused]] void propagate(int gridpoints, double space_width, double start, fftw_complex *psi, double potential(double x), double time_width, int steps
, std::string& output) {
    // create fft and inverse fft plans
    fftw_plan fft = fftw_plan_dft_1d(gridpoints, psi, psi, -1, FFTW_ESTIMATE);
    fftw_plan ifft = fftw_plan_dft_1d(gridpoints, psi, psi, 1, FFTW_ESTIMATE);
    // define potential term (once because it is time-independent)
    fftw_complex V[gridpoints];
    defV(gridpoints, V, potential, space_width, time_width, start);
    // define e^(-i*dt*H) operator (once because it is time-independent);
    fftw_complex H[gridpoints];
    defH(gridpoints, H, space_width, time_width);
    // sign for mult
    int sign;
    // propagation loop
    for (int t = 0; t <= steps; t++) {
        // print completion % to console
        std::cout << GREEN << "\r" << 100*t/steps << "%";
        // print to output file every t_out steps
        if (t % 10 == 0) {
            std::ofstream potwrite;
            potwrite.open(output, std::ios::app);
            if (potwrite.is_open()) {
                for (int i=0; i<gridpoints; i += 128) {
                    potwrite << t * time_width << " " << i*space_width + start << " " << psi[i][0] << " " << psi[i][1] << " " << psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1] << std::endl;
                }
                potwrite.close();
            } else {
                std::cerr << "Failed to open " << output << "." << std::endl;
            }
        }
        // mult space part
        for (int i = 0; i < gridpoints; i++) {
            sign = (i % 2 == 0) ? 1 : -1;
            // sign change for dft
            psi[i][0] *= sign;
            psi[i][1] *= sign;
            double re = psi[i][0];
            double im = psi[i][1];
            // mult by V(x) exponential term
            psi[i][0] = re*V[i][0] - im*V[i][1];
            psi[i][1] = im*V[i][0] + re*V[i][1];
        }
        // execute fft plan
        fftw_execute(fft);
        // multiply psi by e^(-i*dt*H)
        for (int i = 0; i < gridpoints; i++) {
            double re = psi[i][0];
            double im = psi[i][1];
            psi[i][0] = re*H[i][0] - im*H[i][1];
            psi[i][1] = im*H[i][0] + re*H[i][1];
        }
        // execute inverse fft plan
        fftw_execute(ifft);
        // normalize fftw result (fftw uses non-normalized fft algorithm)
        double scale = 1.0 / gridpoints;
        scale_fftw_complex(scale, psi, gridpoints);
        // mult space part
        for (int i = 0; i < gridpoints; i++) {
            sign = (i % 2 == 0) ? 1 : -1;
            // sign change for dft
            psi[i][0] *= sign;
            psi[i][1] *= sign;
            double re = psi[i][0];
            double im = psi[i][1];
            // mult by V(x) exponential term
            psi[i][0] = re*V[i][0] - im*V[i][1];
            psi[i][1] = im*V[i][0] + re*V[i][1];
        }
    }
    std::cout << std::endl;
    // destroy plans
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    for (int i = 0; i < gridpoints; i++) {
        sign = (i % 2 == 0) ? 1 : -1;
        psi[i][0] *= sign;
        psi[i][1] *= sign;
    }
}

// imaginary time propogation
// calculates free-particle operator based on general values and outputs to op
void defiH(int gridpoints, fftw_complex *op, double space_width, double time_width) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = exp(-time_width * psquared(i, gridpoints, space_width) / 2);
        op[i][1] = 0;
    }
}

// creates potential operator array from potential function and outputs to op
void defiV(int gridpoints, fftw_complex *op, double potential(double x), double space_width, double time_width, double start) {
    for (int i = 0; i < gridpoints; i++) {
        op[i][0] = exp(-potential(i * space_width + start) * time_width / 2.0);
        op[i][1] = 0;
    }
}

// propogates wavefunction based on general values
[[maybe_unused]] void ipropagate(int gridpoints, double space_width, double start, fftw_complex *psi, double potential(double x), double time_width, int steps
        , std::string& output) {
    // create fft and inverse fft plans
    fftw_plan fft = fftw_plan_dft_1d(gridpoints, psi, psi, -1, FFTW_ESTIMATE);
    fftw_plan ifft = fftw_plan_dft_1d(gridpoints, psi, psi, 1, FFTW_ESTIMATE);
    // define potential term (once because it is time-independent)
    fftw_complex V[gridpoints];
    defiV(gridpoints, V, potential, space_width, time_width, start);
    // define e^(-i*dt*H) operator (once because it is time-independent);
    fftw_complex H[gridpoints];
    defiH(gridpoints, H, space_width, time_width);
    // sign for mult
    int sign;
    // propagation loop
    for (int t = 0; t <= steps; t++) {
        // print completion % to console
        std::cout << GREEN << "\r" << 100*t/steps << "%";
        // print to output file every t_out steps
        if (t % 10 == 0) {
            std::ofstream potwrite;
            potwrite.open(output, std::ios::app);
            if (potwrite.is_open()) {
                for (int i=0; i<gridpoints; i += 32) {
                    potwrite << t * time_width << " " << i*space_width + start << " " << psi[i][0] << " " << psi[i][1] << " " << psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1] << std::endl;
                }
                potwrite.close();
            } else {
                std::cerr << "Failed to open " << output << "." << std::endl;
            }
        }
        // mult space part
        for (int i = 0; i < gridpoints; i++) {
            sign = (i % 2 == 0) ? 1 : -1;
            // sign change for dft
            psi[i][0] *= sign;
            psi[i][1] *= sign;
            double re = psi[i][0];
            double im = psi[i][1];
            // mult by V(x) exponential term
            psi[i][0] = re*V[i][0] - im*V[i][1];
            psi[i][1] = im*V[i][0] + re*V[i][1];
        }
        // execute fft plan
        fftw_execute(fft);
        // multiply psi by e^(-i*dt*H)
        for (int i = 0; i < gridpoints; i++) {
            double re = psi[i][0];
            double im = psi[i][1];
            psi[i][0] = re*H[i][0] - im*H[i][1];
            psi[i][1] = im*H[i][0] + re*H[i][1];
        }
        // execute inverse fft plan
        fftw_execute(ifft);
        // normalize fftw result (fftw uses non-normalized fft algorithm)
        double scale = 1.0 / gridpoints;
        scale_fftw_complex(scale, psi, gridpoints);
        // mult space part
        for (int i = 0; i < gridpoints; i++) {
            sign = (i % 2 == 0) ? 1 : -1;
            // sign change for dft
            psi[i][0] *= sign;
            psi[i][1] *= sign;
            double re = psi[i][0];
            double im = psi[i][1];
            // mult by V(x) exponential term
            psi[i][0] = re*V[i][0] - im*V[i][1];
            psi[i][1] = im*V[i][0] + re*V[i][1];
        }
    }
    std::cout << std::endl;
    // destroy plans
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    for (int i = 0; i < gridpoints; i++) {
        sign = (i % 2 == 0) ? 1 : -1;
        psi[i][0] *= sign;
        psi[i][1] *= sign;
    }
}

#endif //PROPAGATE_H
