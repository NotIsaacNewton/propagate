//
// Created by Arian Dovald on 6/23/25.
//

#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <functional>

// writes from 1D double function to file
[[maybe_unused]] inline void writeFunction1D(const double& start, const double& end, const double& width,
    const int& gridpoints, const std::string& file, const std::function<double(double)>& function) {
    std::ofstream write;
    write.open(file);
    if (write.is_open()) {
        for (int i=0; i<gridpoints; i++) {
            write << i*width + start << " " << function(i*width + start) << "\n";
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

// reads from file to 1D array
[[maybe_unused]] inline void readArray1D(const std::string& file, std::vector<double>& array) {
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        std::string line;
        int i = 0;
        while (std::getline(read, line)) {
            std::istringstream readline(line);
            double pos; // not used, but needs to moved out of the way
            double val;
            readline >> pos >> val;
            array[i] = val;
            ++i;
        }
        read.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

// writes from 1D array to file
[[maybe_unused]] inline void writeArray1D(const double& start, const double& end, const double& width,
    const std::string& file, const double *function) {
    const double n = (end-start)/width;
    std::ofstream write;
    write.open(file);
    if (write.is_open()) {
        for (int i=0; i<n; i++) {
            write << i*width + start << " " << function[i] << "\n";
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

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
inline inputs readInputs(const std::string& file) {
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        double inputarray[4];
        int intinputarray[4];
        std::print("Reading {}\n",file);
        std::string line;
        int n = 0;
        while (std::getline(read, line)) {
            std::istringstream readline(line);
            (n < 4) ? readline >> inputarray[n] : readline >> intinputarray[n-4];
            ++n;
        }
        read.close();
        // initialize inputs
        inputs in{
            inputarray[0],
            inputarray[1],
            intinputarray[0],
            intinputarray[2],
            inputarray[2],
            inputarray[3],
            intinputarray[1],
            intinputarray[3]
        };
        in.dx = (in.final_pos - in.initial_pos)/in.space_grid;
        in.dt = (in.final_t - in.initial_t)/in.time_grid;
        return in;
    }
    std::cerr << "Failed to open " << file << "." << "\n";
    exit(1);
}

#endif //FILETOOLS_H
