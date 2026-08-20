//
// Created by Arian Dovald on 8/20/26.
//

#include "potentials_td.h"
#include "potentials.h"
#include "console_tools.h"
#include "filetools.h"
#include <iostream>
#include <cmath>
#include <string>

// well + wave
std::function<double(double, double)> electricBarrier(double start_pos, double end_pos,
    double strength_1, double strength_2) {
    return [start_pos, end_pos, strength_1, strength_2](const double x, const double t) {
        return x <= end_pos && x >= start_pos ?
        strength_1 + std::abs(strength_1)*cos(strength_2*x-137/strength_2*t)*exp(-(x-t+4)*(x-t+4)/4) :
        std::abs(strength_1)*cos(strength_2*x-137/strength_2*t)*exp(-(x-t+4)*(x-t+4)/4);
    };
}

// harmonic oscillator + wave
std::function<double(double, double)> electricSho(const double strength, const double k) {
    return [strength, k](const double x, const double t) {
        return strength*x*x/2 + std::abs(strength)*cos(k*x-137/k*t)*exp(-(x-t+4)*(x-t+4)/4);
    };
}

/*std::function<double(double, double)> electricBarrier(double start_pos, double end_pos,
    double strength_1, double strength_2) {
    return [start_pos, end_pos, strength_1, strength_2](const double x, const double t) {
        return x <= end_pos && x >= start_pos ? strength_1 : 0;
    };
}*/

// potential builder + map of options
std::function<double(double, double)> buildPotentialTD(char* argv[]) {
    std::unordered_map<std::string, std::function<std::function<double(double, double)>()>> const potentials = {
        {"electricwell",    [argv] {
            return electricBarrier(std::stod(argv[4]),std::stod(argv[5]),
                std::stod(argv[6]), std::stod(argv[7]));
            }},
        {"electricsho",    [argv] {
            return electricSho(std::stod(argv[6]),std::stod(argv[7]));
            }}
    };
    return potentials.at(argv[3])();
}

// inputs: location/of/input_file location/of/data_directory potential_type pos1 pos2 strength1 strength2
int main(const int argc, char *argv[]) {
    if(argc < 5) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}[location/of/input_file] [location/of/data_directory] [potential_parameters]\n", GREEN);
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
    const double dx = (in.final_pos-in.initial_pos)/(in.space_grid_coarse-1);
    const double dt = (in.final_t-in.initial_t)/(in.time_grid_coarse-1);
    writeFunction2D(in.initial_pos, in.initial_t, dx, dt, in.space_grid_coarse,
        in.time_grid_coarse, potfile, buildPotentialTD(argv));

    // console output
    std::print("{}Potential written!\n\n", GREEN);
    spacerChunky(BLUE);

    return 0;
}
