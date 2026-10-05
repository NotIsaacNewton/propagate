//
// Created by Arian Dovald on 9/29/26.
//

#include "wavepackets_cc.h"
#include "wavepackets.h"

// wavepacket builder for CC
std::vector<std::function<void(double, fftw_complex)>> buildWavepacketCC(const inputs& in) {
    std::vector<std::function<void(double, fftw_complex)>> wavepacket(in.channels);
    wavepacket[0] = buildWavepacket(in);
    for (int c = 1; c < in.channels; c++) {
        wavepacket[c] = zeroState();
    }
    return wavepacket;
}
