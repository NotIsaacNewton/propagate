//
// Created by Arian Dovald on 6/23/25.
//

#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// writes from 1D double function to file
[[maybe_unused]] void writeFunction1D(const double& start, const double& end, const double& width, const std::string& file, double function(double x)) {
    const double n = (end-start)/width;
    std::ofstream write;
    write.open(file);
    if (write.is_open()) {
        for (int i=0; i<n; i++) {
            write << i*width + start << " " << function(i*width + start) << std::endl;
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << std::endl;
    }
}

// writes from 1D array to file
[[maybe_unused]] void writeArray1D(const double& start, const double& end, const double& width, const std::string& file, double *function) {
    const double n = (end-start)/width;
    std::ofstream write;
    write.open(file);
    if (write.is_open()) {
        for (int i=0; i<n; i++) {
            write << i*width + start << " " << function[i] << std::endl;
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << std::endl;
    }
}

void readInputs(double& initial_pos, double& final_pos, int& space_grid,
                double& initial_t, double& final_t, int& time_grid,
                const std::string& file) {
    double inputarray[4];
    int intinputarray[2];
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        std::print("Reading inputs from {}...\n",file);
        std::string line;
        int n = 0;
        while (std::getline(read, line)) {
            std::istringstream readline(line);
            (n < 4) ? readline >> inputarray[n] : readline >> intinputarray[n-4];
            ++n;
        }
        read.close();
        initial_pos = inputarray[0];
        final_pos = inputarray[1];
        space_grid = intinputarray[0];
        initial_t = inputarray[2];
        final_t = inputarray[3];
        time_grid = intinputarray[1];
    } else {
        std::cerr << "Failed to open " << file << "." << std::endl;
    }
}

#endif //FILETOOLS_H
