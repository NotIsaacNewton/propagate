//
// Created by Arian Dovald on 9/2/26.
//

#include <functional>
#include <iostream>
#include <string>
#include <map>
#include <print>
#include <chrono>
#include "file_tools.h"
#include "console_tools.h"
#include "fftw_complex_tools.h"
#include "wavepackets.h"
#include "wavepackets_cc.h"

// inputs: location/of/input_file location/of/data_directory
int main(const int argc, char* argv[]) {
    if (argc != 3) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}[location/of/input_file] [location/of/data_directory]\n", GREEN);
        spacerFancy(RED);
        return 1;
    }

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];

    // introduction
    std::print("\n{}wavepacket\n", BLUE);

    // record the start time
    auto const start = std::chrono::steady_clock::now();

    // spacer
    spacerChunky(BLUE);

    // read input file
    const inputs in = readInputs(inputfile);

    // output psi file
    const std::string psifile = data + "/" + in.input_psi_file;

    // CC wavepacket?
    const bool CC = in.run_type == "TDCC";

    // spacer
    spacer(RESET);

    // write wavefunction
    if (CC) {
        fftw_complex_func_array_to_file(in, data, in.input_psi_file,
            buildWavepacketCC(in));
    } else {
        fftw_complex_func_to_file(in, psifile, buildWavepacket(in));
    }

    // record end time and duration
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    std::print("{}Wavepacket written!\n", GREEN);
    spacerChunky(BLUE);
    std::print("{}Execution time: {} seconds\n\n", BLUE, sec.count());

    return 0;
}
