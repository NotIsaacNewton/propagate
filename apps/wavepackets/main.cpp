//
// Created by Arian Dovald on 9/2/26.
//

#include <functional>
#include <iostream>
#include <string>
#include <map>
#include <print>
#include <chrono>
#include "../../lib/file_tools.h"
#include "../../lib/console_tools.h"
#include "../../lib/fftw_complex_tools.h"
#include "wavepackets.h"

// inputs: location/of/input_file location/of/data_directory wavepacket_type delta momentum position
int main(const int argc, char* argv[]) {
    if (argc < 4) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}[location/of/input_file] [location/of/data_directory] [wavepacket_parameters]\n", GREEN);
        spacerFancy(RED);
        return 1;
    }

    // introduction
    std::print("\n{}wavepacket\n", BLUE);

    // record the start time
    auto const start = std::chrono::steady_clock::now();

    // spacer
    spacerChunky(BLUE);

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string psifile = data+"/psi_initial.dat";

    // read input file
    const inputs in = readInputs(inputfile);

    // spacer
    spacer(RESET);

    // write wavefunction
    fftw_complex_func_to_file(in, psifile, buildWavepacket(argv));

    // record end time and duration
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    std::print("{}Wavepacket written!\n", GREEN);
    spacerChunky(BLUE);
    std::print("{}Execution time: {} seconds\n\n", BLUE, sec.count());

    return 0;
}
