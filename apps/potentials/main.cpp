//
// Created by Arian Dovald on 9/2/26.
//

#include "potentials_td.h"
#include "potentials.h"
#include "../../lib/console_tools.h"
#include "../../lib/file_tools.h"
#include <iostream>
#include <string>
#include <print>
#include <chrono>

// inputs: location/of/input_file location/of/data_directory potential_type pos1 pos2 strength1 strength2 TD/TI
int main(const int argc, char *argv[]) {
    if(argc < 6) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print(
            "{}[location/of/input_file] [location/of/data_directory] [potential_parameters] [TD: TD/TI]\n",
            GREEN);
        spacerFancy(RED);
        return 1;
    }

    // introduction
    std::print("\n{}potential\n", BLUE);

    // record the start time
    auto const start = std::chrono::steady_clock::now();

    // spacer
    spacerChunky(BLUE);

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string potfile = data + "/potential.dat";

    // TD potential?
    const bool potTD = std::string(argv[8]) == "TD";

    // read inputs
    const inputs in = readInputs(inputfile);

    // spacer
    spacer(RESET);

    // write and save potential curve to file
    if (potTD) {
        std::print("Writing potential...\n");
        writeFunction2D(in.initial_pos, in.initial_t, in.dx, in.dt,
            in.space_grid, in.time_grid, potfile, buildPotentialTD(argv));
    } else {
        const double dx = (in.final_pos-in.initial_pos)/(in.space_grid_coarse-1);
        writeFunction1D(in.initial_pos, dx, in.space_grid_coarse,
            potfile, buildPotential(argv));
    }

    // record end time and duration
    auto const end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    spacerChunky("\n" BLUE);
    std::print("{}Execution time: {} seconds\n\n", BLUE, sec.count());

    return 0;
}
