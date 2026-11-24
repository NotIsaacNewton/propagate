//
// Created by Arian Dovald on 9/28/26.
//

#ifndef PROPAGATE_POTENTIALS_TD_CC_H
#define PROPAGATE_POTENTIALS_TD_CC_H

#include <functional>
#include "matrix_tools.h"
#include "file_tools.h"
#include "potentials.h"
#include <complex>

// harmonic oscillator coupled to decaying exponential by constant
double SHO(double x);
double expDecay(double x, int c, int channels, double range);
std::function<hermitian_matrix(double, double)> coupledSHO(const inputs& in);

// this is called for the potentials during hybrid runs
std::function<std::vector<double>(double, int)> hybridPotentialFunctions(int states, double dE);
std::function<hermitian_matrix(double, double)> hybridPotential(const inputs& in);

// potential from file: assumes one bound ground state
// requires a file containing the bound (ground state) potential curve
// requires N files containing potentials curves for N residual states
// requires N files containing bound-continuum dipole couplings for all N residual states, at M continuum energies
// requires N*(N-1)/2 files with M^2 continuum-continuum dipole couplings (since each potential comes with its
//  own continuum) ordered by rows of increasing energy (E_00 -> E_0M, E_10 -> E_1M, ..., E_M0 -> E_MM)
// grabs all bound + residual potential curves and interpolates each one
struct potentials_from_file {
    double dx{}; // width of coarse potential grid
    std::vector<double> grid; // stores grid on which potential is defined
    std::vector<std::vector<double>> potentials; // coarse potential array
    std::vector<spline_interp> interpolators;
    // constructor
    explicit potentials_from_file(const inputs& in) : dx((in.final_pos - in.initial_pos) / (in.space_grid_coarse - 1)),
                                                   grid(in.space_grid_coarse),
                                                   potentials(in.bound_channels+in.residual_channels,
                                                       std::vector<double>(in.space_grid_coarse)),
                                                   interpolators([&] {
                                                       buildGrid(in, grid, dx);
                                                       getPotentials(in);
                                                       return interpPotentials(in, grid, potentials);
                                                   }()) {}
    void getPotentials(const inputs& in); // gets potentials
    static std::vector<spline_interp> interpPotentials(const inputs& in, const std::vector<double>& grid,
    const std::vector<std::vector<double>>& potentials); // interpolates potentials
    // accessor (read-only)
    double operator()(const int c, const double x) const {
        return interpolators[c].interp_threaded(x);
    }
    // no copying or moving potentialsFromFile objects: spline_interp spans would not follow along
    potentials_from_file(const potentials_from_file&) = delete;
    potentials_from_file& operator=(const potentials_from_file&) = delete;
    potentials_from_file(potentials_from_file&&) = delete;
    potentials_from_file& operator=(potentials_from_file&&) = delete;
};
// grabs all couplings
struct couplings_from_file {
    hermitian_matrix coup; // holds all coupling operators
    // constructor
    explicit couplings_from_file(const inputs& in) :
    coup(in.channels) {
        couplingsFromFileBtoC(in);
        couplingsFromFileCtoC(in);
    }
    // grabs all b-c dipole couplings into hermitian_matrix coup
    void couplingsFromFileBtoC(const inputs& in);
    // grabs all c-c dipole couplings into hermitian_matrix coup
    void couplingsFromFileCtoC(const inputs& in);
    // accessor (read-only)
    const std::complex<double>& operator()(const int i, const int j) const { return coup(i, j); }
};
// uses the above methods to get a complete set of interpolators for all TDCC objects that can be packaged
//  into the standard std::function<hermitian_matrix(double, double)>
struct potential_from_file_tdcc {
    const potentials_from_file pots; // holds all potentials
    const couplings_from_file coups; // holds all couplings
    const inputs in;
    explicit potential_from_file_tdcc(const inputs& in) : pots(in), coups(in), in(in) {}
    // accessor (read-only)
    hermitian_matrix operator()(const double x, const double) const {
        hermitian_matrix pot_matrix(in.channels);
        for (int c = 0; c < in.bound_channels; c++) {
            pot_matrix(c, c) = std::complex(pots(c, x), 0.0);
        }
        for (int c = in.bound_channels; c < in.bound_channels + in.residual_channels; c++) {
            for (int cont = 0; cont < in.continuum_channels; cont++) {
                pot_matrix(in.bound_channels + (c - in.bound_channels)*in.continuum_channels + cont,
                    in.bound_channels + (c - in.bound_channels)*in.continuum_channels + cont) =
                        std::complex(pots(c, x) + cont*in.dE, 0.0);
            }
        }
        for (int c1 = 0; c1 < in.channels; c1++) {
            for (int c2 = c1+1; c2 < in.channels; c2++) {
                pot_matrix(c1,c2) = coups(c1, c2);
            }
        }
        return pot_matrix;
    }
};
// uses the above methods to get a complete set of interpolators for all TDCC objects that can be packaged
//  into the standard std::function<hermitian_matrix(double, double)>
struct hybrid_from_file_tdcc {
    const inputs in;
    const std::function<std::vector<double>(double, int)> potentials;
    const couplings_from_file coups; // holds all couplings
    explicit hybrid_from_file_tdcc(const inputs& in) : in(in),
    potentials(hybridPotentialFunctions(in.bound_channels + in.residual_channels, in.dE)), coups(in) {}
    // accessor (read-only)
    hermitian_matrix operator()(const double x, const double) const {
        hermitian_matrix pot_matrix(in.channels);
        for (int c = 0; c < in.bound_channels; c++) {
            pot_matrix(c, c) = potentials(x, 0)[c];
        }
        for (int c = in.bound_channels; c < in.bound_channels + in.residual_channels; c++) {
            for (int cont = 0; cont < in.continuum_channels; cont++) {
                pot_matrix(in.bound_channels + (c - in.bound_channels)*in.continuum_channels + cont,
                    in.bound_channels + (c - in.bound_channels)*in.continuum_channels + cont) =
                        potentials(x, cont)[c];
            }
        }
        for (int c1 = 0; c1 < in.channels; c1++) {
            for (int c2 = c1+1; c2 < in.channels; c2++) {
                pot_matrix(c1,c2) = coups(c1, c2);
            }
        }
        return pot_matrix;
    }
};

// potential builder + map of options
std::function<hermitian_matrix(double, double)> buildPotentialTDCC(const inputs& in);

#endif //PROPAGATE_POTENTIALS_TD_CC_H
