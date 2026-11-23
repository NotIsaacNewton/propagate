//
// Created by Arian Dovald on 9/8/25.
//

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <functional>
#include "file_tools.h"
#include "interpolate_1d.h"

// potential functions
// step potential
std::function<double(double)> step(double pos, double strength_1, double strength_2);
// harmonic oscillator
std::function<double(double)> sho(double strength);
// well
std::function<double(double)> barrier(double start_pos, double end_pos, double strength);
// triangle
std::function<double(double)> triangle(double strength, double pos);
// wall
std::function<double(double)> wall(double pos, double strength);
// soft-core hydrogen
std::function<double(double)> hydrogen();

// builds a coarse grid
void buildGrid(const inputs& in, std::vector<double>& grid, const double& dx);

// potential from filer
// reads in potential and holds interpolation object with associated data
struct potential_from_file {
    double dx{}; // width of coarse potential grid
    std::vector<double> grid; // stores grid on which potential is defined
    std::vector<double> potential; // coarse potential array
    spline_interp interpolator;
    // constructor
    explicit potential_from_file(const inputs& in) : dx((in.final_pos - in.initial_pos) / (in.space_grid_coarse - 1)),
                                                   grid(in.space_grid_coarse),
                                                   potential(in.space_grid_coarse),
                                                   interpolator([&] {
                                                       buildGrid(in, grid, dx);
                                                       readArray1D(in.pot_read_from_file, potential);
                                                       return spline_interp(grid, potential);
                                                   }()) {}
    // accessor (read-only)
    double operator()(const double x) const {
        return interpolator.interp_threaded(x);
    }
    // no copying or moving potentialFromFile objects: spline_interp spans would not follow along
    potential_from_file(const potential_from_file&) = delete;
    potential_from_file& operator=(const potential_from_file&) = delete;
    potential_from_file(potential_from_file&&) = delete;
    potential_from_file& operator=(potential_from_file&&) = delete;
};

// map of options
std::function<double(double)> buildPotential(const inputs& in);

#endif //POTENTIALS_H