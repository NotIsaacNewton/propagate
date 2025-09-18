//
// Created by Arian Dovald on 6/30/25.
//

#ifndef WAVEPACKETS_H
#define WAVEPACKETS_H

#include <unordered_map>
#include <functional>
#include <string>
#include "fftw3.h"

// Gaussian wave-packet
// Gaussian
double gaussian(double x, double delta, double pos);
// gives Gaussian momentum
void gaussianMoving(double x, double delta, double momentum, double pos, fftw_complex out);
// prepares wave-packet
std::function<void(double, fftw_complex)> gaussianWP(double delta, double momentum, double pos);

// SHO ground-state
std::function<void(double, fftw_complex)> shoGround();

// SHO first excited state
std::function<void(double, fftw_complex)> shoExcited();

// test state (for imaginary time propagation)
std::function<void(double, fftw_complex)> test();

// map of options
std::unordered_map<std::string, std::function<void(double, fftw_complex)>> wpOptions(char* argv[]);

#endif //WAVEPACKETS_H
