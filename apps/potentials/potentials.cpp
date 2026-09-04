//
// Created by Arian Dovald on 9/9/25.
//

#include <unordered_map>
#include "potentials.h"
#include "file_tools.h"

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
std::function<double(double)> triangle(const double strength, const double pos) {
    return [strength, pos](const double x) {
        return x > -pos && x <= 0 ? -strength*x-strength : x > 0 && x <= pos ? strength*x-strength : 0;
    };
}
// wall
std::function<double(double)> wall(const double pos, const double strength) {
    return [pos, strength](const double x) {
        return x <= pos ? 0 : strength;
    };
}

// potential builder + map of options
std::function<double(double)> buildPotential(char* argv[]) {
    std::unordered_map<std::string, std::function<std::function<double(double)>()>> const potentials = {
        {"step",    [argv] {
            return step(std::stod(argv[4]),std::stod(argv[6]),std::stod(argv[7]));
            }},
        {"sho",     [argv] {
            return sho(std::stod(argv[6]));
            }},
        {"well",    [argv] {
            return barrier(std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]));
            }},
        {"triangle",[argv] {
            return triangle(std::stod(argv[6]),std::stod(argv[4]));
            }},
        {"wall",    [argv] {
            return wall(std::stod(argv[4]),std::stod(argv[6]));
            }}
    };
    return potentials.at(argv[3])();
}
