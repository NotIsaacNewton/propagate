//
// Created by Arian Dovald on 9/8/25.
//

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <unordered_map>
#include <functional>

// potential functions
// step potential
std::function<double(double)> step(double pos, double strength_1, double strength_2);
// harmonic oscillator
std::function<double(double)> sho(double strength);
// well
std::function<double(double)> barrier(double start_pos, double end_pos, double strength);
// triangle
std::function<double(double)> triangle(double strength);
// wall
std::function<double(double)> wall(double pos, double strength);

// map of options
std::unordered_map<std::string, std::function<double(double)>> potOptions(char* argv[]);

#endif //POTENTIALS_H