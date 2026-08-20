//
// Created by Arian Dovald on 9/18/25.
//

#include "filetools.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <print>

// writes from 1D double function to file
void writeFunction1D(const double& start, const double& width,
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
void readArray1D(const std::string& file, std::vector<double>& array) {
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
void writeArray1D(const double& start, const double& end, const double& width, const int& gridpoints,
    const std::string& file, const std::vector<double>& function) {
    std::ofstream write;
    write.open(file);
    if (write.is_open()) {
        for (int i=0; i<gridpoints; i++) {
            write << i*width + start << " " << function[i] << "\n";
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

// writes from 2D double function to file
void writeFunction2D(const double& start_x, const double& start_y, const double& dx, const double& dy,
    const int& width, const int& height, const std::string& file,
    const std::function<double(double, double)>& function) {
    std::ofstream write;
    write.open(file);
    if (write.is_open()) {
        for (int i=0; i<height; i++) {
            for (int j=0; j<width; j++) {
                write << function(j*dx + start_x, i*dy + start_y) << " ";
            }
            write << "\n";
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

// reads from file to 2D array
void readArray2D(const std::string& file, std::vector<std::vector<double>>& array,
    const int& width, const int& height) {
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        array.assign(height, std::vector<double>(width));
        std::string line;
        for (int i=0; i<height; i++) {
            std::getline(read, line);
            std::istringstream readline(line);
            for (int j=0; j<width; j++) {
                readline >> array[i][j];
            }
        }
        read.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

// writes from 2D array to file
void writeArray2D(const std::string& file, const std::vector<std::vector<double>>& array,
    const int& width, const int& height) {
    std::ofstream write;
    write.open(file);
    if (write.is_open()) {
        for (int i=0; i<height; i++) {
            for (int j=0; j<width; j++) {
                write << array[i][j] << " ";
            }
            write << "\n";
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

// read inputs from file
inputs readInputs(const std::string& file) {
    std::ifstream read;
    read.open(file);
    if (read.is_open()) {
        double inputarray[4];
        int intinputarray[6];
        std::print("Reading {}\n",file);
        std::string line;
        int n = 0;
        while (std::getline(read, line)) {
            std::istringstream readline(line);
            n < 4 ? readline >> inputarray[n] : readline >> intinputarray[n-4];
            ++n;
        }
        read.close();
        // initialize inputs
        inputs in{
            inputarray[0],
            inputarray[1],
            intinputarray[0],
            intinputarray[3],
            intinputarray[1],
            inputarray[2],
            inputarray[3],
            intinputarray[2],
            intinputarray[4],
            intinputarray[5]
        };
        in.dx = (in.final_pos - in.initial_pos)/(in.space_grid - 1);
        in.dt = (in.final_t - in.initial_t)/(in.time_grid-1);
        return in;
    }
    std::cerr << "Failed to open " << file << "." << "\n";
    exit(1);
}
