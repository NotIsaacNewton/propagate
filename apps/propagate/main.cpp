//
// Created by Arian Dovald on 9/2/26.
//

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <print>
#include "../../lib/file_tools.h"
#include "../../lib/console_tools.h"
#include "propagate.h"
#include "propagate_td.h"

// inputs: location/of/input_file location/of/data_directory TD/improp
int main(const int argc, const char* argv[]) {
    if (argc != 4) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print(
            "{}[location/of/input_file] [location/of/data_directory] [TD/true/false] \n",
            GREEN
            );
        spacerFancy(RED);
        return 1;
    }

    // introduction
    std::print("\n{}propagate\n", BLUE);

    // record the start time
    auto start = std::chrono::steady_clock::now();

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string psiout = data + "/psi_final.dat";

    // imaginary propagation?
    const bool imProp = std::string(argv[3]) == "true";
    // TD propagation?
    const bool propTD = std::string(argv[3]) == "TD";

    // spacer
    spacerChunky(BLUE);

    // read input file
    const inputs in = readInputs(inputfile);

    // warning if odd number of gridpoints
    if (in.space_grid % 2 != 0) {
        spacerFancy(YELLOW);
        std::cerr << YELLOW << "WARNING: Odd number of gridpoints may produce asymmetric momentum grid.\n" << RESET;
        spacerFancy(YELLOW);
    }
    // warning if large dt
    if (in.time_grid < 1000) {
        spacerFancy(YELLOW);
        std::cerr << YELLOW << "WARNING: Large dt may fail to resolve rapid dynamics. Use dt < pi/E_max.\n" << RESET;
        spacerFancy(YELLOW);
    }

    // clear output file
    std::ofstream clearout;
    clearout.open(psiout);
    clearout.close();

    // spacer
    spacerThick(RESET);

    // propagate wf
    if (propTD) { propagateTD(in, data); } else { propagate(in, data, imProp); }

    // spacer
    spacerThick(RESET);

    // record end time and duration
    auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    std::print("Done :)\n");
    spacerChunky(BLUE);
    std::print("{}Execution time: {} seconds\n\n", BLUE, sec.count());

    return 0;
}
