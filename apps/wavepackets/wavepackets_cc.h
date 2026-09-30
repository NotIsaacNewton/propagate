//
// Created by Arian Dovald on 9/29/26.
//

#ifndef PROPAGATE_WAVEPACKETS_CC_H
#define PROPAGATE_WAVEPACKETS_CC_H

#include <functional>
#include "fftw3.h"



// wavepacket builder for CC
std::vector<std::function<void(double, fftw_complex)>> buildWavepacketCC(char* argv[], int channels);

#endif //PROPAGATE_WAVEPACKETS_CC_H
