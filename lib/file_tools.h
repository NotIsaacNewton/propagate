//
// Created by Arian Dovald on 6/23/25.
//

#ifndef FILE_TOOLS_H
#define FILE_TOOLS_H

#include <fstream>
#include <string>
#include <functional>
#include "matrix_tools.h"

// writes from 1D double function to file
void writeFunction1D(const double& start, const double& width,
    const int& gridpoints, const std::string& file, const std::function<double(double)>& function);

// reads from file to 1D array
void readArray1D(const std::string& file, std::vector<double>& array);

// writes from 1D array to file
void writeArray1D(const double& start, const double& width, const int& gridpoints,
    const std::string& file, const std::vector<double>& function);

// writes from 2D double function to file
void writeFunction2D(const double& start_x, const double& start_y, const double& dx, const double& dy,
    const int& width, const int& height, const std::string& file,
    const std::function<double(double, double)>& function);

// writes from 2D double function of hermitian matrices to file
void writeHermitian2D(const double& start_x, const double& start_y, const double& dx, const double& dy,
    const int& width, const int& height, const std::string& file,
    const std::function<hermitian_matrix(double, double)>& function);

// reads from file to 2D array
void readArray2D(const std::string& file, std::vector<std::vector<double>>& array,
    const int& width, const int& height);

// writes from 2D array to file
void writeArray2D(const std::string& file, const std::vector<std::vector<double>>& array,
    const int& width, const int& height);

// TODO: make readInputs define all new fields here
//  modify calculate.sh to write inputs.txt with all fields below
//  document all input options and behavior
// inputs go here
struct inputs {
    // main file names
    std::string input_psi_file;
    std::string output_psi_file;
    std::string pot_file;
    // system data
    std::string run_type;
    int channels;
    // space data
    double initial_pos;
    double final_pos;
    int space_grid;
    int nx_prints;
    int space_grid_coarse;
    double dx;
    // time data
    double initial_t;
    double final_t;
    int time_grid;
    int nt_prints;
    double dt;
    // potential data
    std::string potential_type;
    std::string pot_read_from_file;
    double pot_position_1;
    double pot_position_2;
    double strength_1;
    double strength_2;
    // wavepacket data
    std::string wavepacket_type;
    std::string wp_read_from_file;
    double delta;
    double momentum;
    double wp_position;
};

// read inputs from file
inputs readInputs(const std::string& file);

// wavefunction + buffer struct
struct wfOutput {
    std::ofstream wf;
    std::vector<double> buffer;
};

// opens wavefunction output file
wfOutput openWFOutputFile(const inputs& in, const std::string &data);

// opens wavefunction output file, but with a larger buffer (by a factor of # of channels)
wfOutput openWFOutputFileCC(const inputs& in, const std::string &data);

#endif //FILE_TOOLS_H
