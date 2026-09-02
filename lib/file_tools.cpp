//
// Created by Arian Dovald on 9/18/25.
//

#include "file_tools.h"
#include "console_tools.h"
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
void writeArray1D(const double& start, const double& width, const int& gridpoints,
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
    std::ofstream write(file, std::ios::binary);
    if (!write.is_open()) {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
    if (write.is_open()) {
        std::vector<double> temp(width);
        for (int i=0; i<height; i++) {
            !((i+1) % (height / 10)) ? progressBar(GREEN, 100*(i+1)/height) : reset();
            for (int j=0; j<width; j++) {
                temp[j] = function(j*dx + start_x, i*dy + start_y);
            }
            write.write(
            reinterpret_cast<const char*>(temp.data()), static_cast<std::streamsize>(width * sizeof(double)));
        }
        write.close();
    } else {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
}

// reads from file to 2D array
void readArray2D(const std::string& file, std::vector<std::vector<double>>& array,
    const int& width, const int& height) {
    if (std::ifstream read(file, std::ios::binary); read.is_open()) {
        array.assign(height, std::vector<double>(width));
        for (int i=0; i<height; i++) {
            read.read(reinterpret_cast<char*>(array[i].data()), static_cast<std::streamsize>(width * sizeof(double)));
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
            .initial_pos = inputarray[0],
            .final_pos = inputarray[1],
            .space_grid = intinputarray[0],
            .nx_prints = intinputarray[3],
            .space_grid_coarse = intinputarray[1],
            .initial_t = inputarray[2],
            .final_t = inputarray[3],
            .time_grid = intinputarray[2],
            .nt_prints = intinputarray[4],
            .time_grid_coarse = intinputarray[5]
        };
        in.dx = (in.final_pos - in.initial_pos)/(in.space_grid - 1);
        in.dt = (in.final_t - in.initial_t)/(in.time_grid-1);
        std::print("Read {}\n",file);
        return in;
    }
    std::cerr << "Failed to open " << file << "." << "\n";
    exit(1);
}

// opens wavefunction output file
wfOutput openWFOutputFile(const inputs& in, const std::string &data) {
    const std::string output = data + "/psi_final.dat";
    std::ofstream wf(output, std::ios::app | std::ios::binary);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // prepare output buffer for entire set of points
    std::vector<double> buffer;
    buffer.reserve((in.time_grid / in.nt_prints + 1) * (in.space_grid / in.nx_prints) * 2);
    return {.wf = std::move(wf), .buffer = buffer};
}
