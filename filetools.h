//
// Created by Arian Dovald on 6/23/25.
//

#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <string>
#include <functional>

// writes from 1D double function to file
void writeFunction1D(const double& start, const double& width,
    const int& gridpoints, const std::string& file, const std::function<double(double)>& function);

// reads from file to 1D array
void readArray1D(const std::string& file, std::vector<double>& array);

// writes from 1D array to file
void writeArray1D(const double& start, const double& end, const double& width,
    const std::string& file, const double *function);

// inputs go here
struct inputs {
    // space stuff
    double initial_pos;
    double final_pos;
    int space_grid;
    int nx_prints;
    // time stuff
    double initial_t;
    double final_t;
    int time_grid;
    int nt_prints;
    // space-time widths
    double dx;
    double dt;
};

// read inputs from file
inputs readInputs(const std::string& file);

#endif //FILETOOLS_H
