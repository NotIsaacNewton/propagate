//
// Created by Arian Dovald on 9/9/25.
//

#include <unordered_map>
#include "potentials.h"

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
// soft-core hydrogen
std::function<double(double)> hydrogen() {
    return [=](const double x) {
        return -1/std::sqrt(1 + x*x);
    };
}

// builds a coarse grid
void buildGrid(const inputs& in, std::vector<double>& grid, const double& dx) {
    // write coarse potential grid
    #pragma omp parallel for default(none) shared(in, grid, dx)
    for (int i = 0; i < in.space_grid_coarse; i++) {
        grid[i] = in.initial_pos + i * dx;
    }
}

// potential builder + map of options
std::function<double(double)> buildPotential(const inputs& in) {
    std::unordered_map<std::string, std::function<std::function<double(double)>()>> const potentials = {
        {"step",    [in] {
            return step(in.pot_position_1,in.strength_1,in.strength_2);
            }},
        {"sho",     [in] {
            return sho(in.strength_1);
            }},
        {"well",    [in] {
            return barrier(in.pot_position_1,in.pot_position_2,in.strength_1);
            }},
        {"triangle",[in] {
            return triangle(in.strength_1,in.pot_position_1);
            }},
        {"wall",    [in] {
            return wall(in.pot_position_1,in.strength_1);
            }},
        {"hydrogen",  [] {
            return hydrogen();
            }},
        {"file",    [in] {
            auto pot = std::make_shared<potential_from_file>(in);
            return [pot](const double x) {
                return pot->operator()(x);
            };
            }}
    };
    return potentials.at(in.potential_type)();
}
