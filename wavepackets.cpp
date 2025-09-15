//
// Created by Arian Dovald on 9/9/25.
//

#include <functional>
#include <map>
#include "filetools.h"
#include "fftw_complex_tools.h"
#include "wavepackets.h"
#include "console_tools.h"

// inputs: location/of/input_file location/of/data_directory wavepacket_type delta momentum position
int main(const int argc, char* argv[]) {
    if (argc < 4) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}location/of/input_file location/of/data_directory wavepacket_parameters\n", GREEN);
        spacerFancy(RED);
        return 1;
    }

    // spacer
    spacerChunky(BLUE);
    std::print("\n");

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string psifile = data+"/psi_initial.dat";

    // read input file
    const inputs in = readInputs(inputfile);

    // spacer
    spacer(RESET);

    // write wavefunction
    fftw_complex_func_to_file(in, psifile, wpOptions(argv).at(argv[3]));

    // console output
    std::print("{}Wavepacket written!\n\n", GREEN);
    spacerChunky(BLUE);

    return 0;
}
