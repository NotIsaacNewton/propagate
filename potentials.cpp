//
// Created by Arian Dovald on 9/9/25.
//

#include "potentials.h"
#include "fftw_complex_tools.h"
#include "filetools.h"
#include "console_tools.h"

// TODO: add time-dependent potentials

// inputs: location/of/input_file location/of/data_directory potential_type pos1 pos2 strength1 strength2
int main(const int argc, char *argv[]) {
    if(argc < 5) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}location/of/input_file location/of/data_directory potential_parameters\n", GREEN);
        spacerFancy(RED);
        return 1;
    }

    // spacer
    spacerChunky(BLUE);
    std::print("\n");

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string potfile = data + "/potential.dat";

    // read inputs
    const inputs in = readInputs(inputfile);

    // spacer
    spacer(RESET);

    // write and save potential curve to file

    writeFunction1D(in.initial_pos, in.final_pos, in.dx, in.space_grid,
        potfile, potOptions(argv).at(argv[3]));

    // console output
    std::print("{}Potential written!\n\n", GREEN);
    spacerChunky(BLUE);

    return 0;
}
