//
// Created by Arian Dovald on 9/9/25.
//

#include <iostream>
#include "potentials.h"
#include "fftw_complex_tools.h"
#include "filetools.h"
#include "console_tools.h"

// TODO: add time-dependent potentials

// potential functions
// step potential
std::function<double(double)> step(const double pos,
    const double strength_1, const double strength_2) {
    return [pos, strength_1, strength_2](const double x) {
        return x <= pos ? strength_1 : strength_2;
    };
}
// harmonic oscillator
std::function<double(double)> sho(const double strength) {
    return [strength](const double x) {
        return strength*x*x/2;
    };
}
// well
std::function<double(double)> barrier(const double start_pos,
    const double end_pos, const double strength) {
    return [start_pos, end_pos, strength](const double x) {
        return x <= end_pos && x >= start_pos ? strength : 0;
    };
}
// triangle
std::function<double(double)> triangle(const double strength) {
    return [strength](const double x) {
        return x > -3 && x <= 0 ? -strength*x : x > 0 && x <= 3 ? strength*x : 0;
    };
}
// wall
std::function<double(double)> wall(const double pos, const double strength) {
    return [pos, strength](const double x) {
        return x <= pos ? 0 : strength;
    };
}

// map of options
std::unordered_map<std::string, std::function<double(double)>> potOptions(char* argv[]) {
    return {
            {
                {"step", step(std::stod(argv[4]),std::stod(argv[6]),std::stod(argv[7]))},
                {"sho", sho(std::stod(argv[6]))},
                {"well", barrier(std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]))},
                {"triangle", triangle(std::stod(argv[6]))},
                {"wall", wall(std::stod(argv[4]),std::stod(argv[6]))}
            }
    };
}

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

    writeFunction1D(in.initial_pos, in.dx, in.space_grid,
        potfile, potOptions(argv).at(argv[3]));

    // console output
    std::print("{}Potential written!\n\n", GREEN);
    spacerChunky(BLUE);

    return 0;
}
