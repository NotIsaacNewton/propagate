//
// Created by Arian Dovald on 9/2/26.
//

#include "potentials_td_cc.h"
#include "potentials_td.h"
#include "potentials.h"
#include "console_tools.h"
#include "file_tools.h"
#include <iostream>
#include <string>
#include <print>
#include <chrono>

// inputs: location/of/input_file location/of/data_directory
int main(const int argc, char *argv[]) {
    if(argc != 3) {
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
    std::print("\n{}potential\n", BLUE);

    // record the start time
    auto const start = std::chrono::steady_clock::now();

    // spacer
    spacerChunky(BLUE);

    // read inputs
    const inputs in = readInputs(inputfile);

    // potential file
    const std::string potfile = data + "/" + in.pot_file;

    // TD potential?
    const bool potTD = in.run_type == "TD";
    // TDCC potential?
    const bool potTDCC = in.run_type == "TDCC";

    // spacer
    spacer(RESET);

    // write and save potential curve to file
    if (potTD) {
        std::print("Writing potential...\n");
        writeFunction2D(in.initial_pos, in.initial_t, in.dx, in.dt,
            in.space_grid, in.time_grid, potfile, buildPotentialTD(in));
    } else if (potTDCC) {
        std::print("Writing potential...\n");
        writeHermitian2D(in.initial_pos, in.initial_t, in.dx, in.dt,
            in.space_grid, in.time_grid, potfile, buildPotentialTDCC(in));
    } else {
        const double dx = (in.final_pos-in.initial_pos)/(in.space_grid_coarse-1);
        writeFunction1D(in.initial_pos, dx, in.space_grid_coarse,
            potfile, buildPotential(in));
    }

    // record end time and duration
    auto const end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    spacerChunky("\n" BLUE);
    std::print("{}Execution time: {} seconds\n\n", BLUE, sec.count());

    return 0;
}
