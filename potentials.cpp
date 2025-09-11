//
// Created by Arian Dovald on 9/9/25.
//

#include "potentials.h"
#include "fftw_complex_tools.h"
#include "filetools.h"
#include "console_tools.h"

// TODO: generalize inputs for a wider range of potentials

// inputs: location/of/input_file location/of/data_directory
int main(int argc, char *argv[]) {
    if(argc != 3) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}location/of/input_file location/of/data_directory\n", GREEN);
        spacerFancy(RED);
    }

    // spacer
    spacerChunky(BLUE);
    std::print("\n");

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string potfile = data + "/potential.dat";

    // read inputs
    inputs in = readInputs(inputfile);

    // spacer
    spacer(RESET);

    // temporary fftw_complex array
    const auto temp = fftw_alloc_complex(in.space_grid);
    std::unique_ptr<fftw_complex, void(*)(void*)> psip{temp, fftw_free};
    // write and save potential operator to file
    definePotentialOperator(in.space_grid, temp, sho(50),
        in.dx, in.dt, in.initial_pos);
    fftw_complex_array_to_file(in.initial_pos, in.final_pos, in.dx, potfile, temp);

    // console output
    std::print("{}Potential written!\n\n", GREEN);
    spacerChunky(BLUE);

    return 1;
}
