//
// Created by Arian Dovald on 6/30/25.
//

#ifndef WAVEPACKETS_H
#define WAVEPACKETS_H

#include <functional>
#include "file_tools.h"
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

// utility function for filling intially-empty channels
std::function<void(double, fftw_complex)> zeroState();

// read function from file of raw doubles
std::function<void(double, fftw_complex)> psiFromFile(const inputs& in);

// map of options
std::function<void(double, fftw_complex)> buildWavepacket(const inputs& in);

#endif //WAVEPACKETS_H
