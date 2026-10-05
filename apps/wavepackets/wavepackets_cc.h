//
// Created by Arian Dovald on 9/29/26.
//

#ifndef PROPAGATE_WAVEPACKETS_CC_H
#define PROPAGATE_WAVEPACKETS_CC_H

#include <functional>
#include "fftw3.h"
#include "file_tools.h"

// wavepacket builder for CC
std::vector<std::function<void(double, fftw_complex)>> buildWavepacketCC(char* argv[], const inputs& in, int channels);

#endif //PROPAGATE_WAVEPACKETS_CC_H
