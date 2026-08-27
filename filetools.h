//
// Created by Arian Dovald on 6/23/25.
//

#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <fstream>
#include <string>
#include <functional>

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

// reads from file to 2D array
void readArray2D(const std::string& file, std::vector<std::vector<double>>& array,
    const int& width, const int& height);

// writes from 2D array to file
void writeArray2D(const std::string& file, const std::vector<std::vector<double>>& array,
    const int& width, const int& height);

// inputs go here
struct inputs {
    // space stuff
    double initial_pos;
    double final_pos;
    int space_grid;
    int nx_prints;
    int space_grid_coarse;
    // time stuff
    double initial_t;
    double final_t;
    int time_grid;
    int nt_prints;
    int time_grid_coarse;
    // space-time widths
    double dx;
    double dt;
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

#endif //FILETOOLS_H
