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
        strength_1 + std::abs(strength_1)*cos(strength_2*x-137*strength_2*t)*exp(-(x-t+8)*(x-t+8)/16) :
        std::abs(strength_1)*cos(strength_2*x-137*strength_2*t)*exp(-(x-t+8)*(x-t+8)/16);
    };
}

// harmonic oscillator + wave
std::function<double(double, double)> electricSho(const double strength, const double k) {
    return [strength, k](const double x, const double t) {
        return strength*x*x/2 + std::abs(strength)*50*cos(k*x-137*k*t)*exp(-(x-t+4)*(x-t+4)/1);
    };
}

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

    // introduction
    std::print("\n{}potential_td\n", BLUE);

    // record the start time
    auto const start = std::chrono::steady_clock::now();

    // spacer
    spacerChunky(BLUE);

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string potfile = data + "/potential.dat";

    // read inputs
    const inputs in = readInputs(inputfile);

    // spacer
    spacer(RESET);

    // write and save potential curve to file
    writeFunction2D(in.initial_pos, in.initial_t, in.dx, in.dt,
        in.space_grid, in.time_grid, potfile, buildPotentialTD(argv));

    // record end time and duration
    auto const end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    std::print("{}Potential written!\n", GREEN);
    spacerChunky(BLUE);
    std::print("{}Execution time: {} seconds\n\n", BLUE, sec.count());

    return 0;
}
