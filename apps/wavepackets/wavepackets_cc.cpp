//
// Created by Arian Dovald on 9/29/26.
//

#include "wavepackets_cc.h"
#include "wavepackets.h"

// wavepacket builder for CC
std::vector<std::function<void(double, fftw_complex)>> buildWavepacketCC(char* argv[], const inputs& in,
    const int channels) {
    std::vector<std::function<void(double, fftw_complex)>> wavepacket(channels);
    wavepacket[0] = buildWavepacket(argv, in);
    for (int c = 1; c < channels; c++) {
        wavepacket[c] = zeroState();
    }
    return wavepacket;
}
