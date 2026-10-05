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
            #pragma omp parallel for default(none) shared(i, temp, width, start_x, start_y, dx, dy, function)
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

// writes from 2D hermitian matrix to file
void writeHermitian2D(const double& start_x, const double& start_y, const double& dx, const double& dy,
    const int& width, const int& height, const std::string& file,
    const std::function<hermitian_matrix(double, double)>& function) {
    std::ofstream write(file, std::ios::binary);
    if (!write.is_open()) {
        std::cerr << "Failed to open " << file << "." << "\n";
    }
    if (write.is_open()) {
        const int channels = function(start_x, start_y).N;
        std::vector<double> temp(width*(channels*channels+channels));
        for (int i = 0; i < height; i++) {
            !((i+1) % (height / 10)) ? progressBar(GREEN, 100*(i+1)/height) : reset();
            #pragma omp parallel for default(none) shared(i, width, channels, function, start_x, start_y, dx, dy, temp)
            for (int j = 0; j < width; j++) {
                const hermitian_matrix matrix = function(j*dx + start_x, i*dy + start_y);
                for (int c1 = 0; c1 < channels; c1++) {
                    for (int c2 = c1; c2 < channels; c2++) {
                        const auto value = matrix(c1,c2);
                        temp[2*(j + c1*(channels-(c1+1)/2)*width + c2*width)]   = value.real();
                        temp[2*(j + c1*(channels-(c1+1)/2)*width + c2*width)+1] = value.imag();
                    }
                }
            }
            write.write(
                reinterpret_cast<const char*>(temp.data()),
                static_cast<std::streamsize>(temp.size() * sizeof(double)));
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
            !((i+1) % (height / 10)) ? progressBar(GREEN, 100*(i+1)/height) : reset();
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
        std::string strings[8];
        double doubles[11];
        int ints[6];
        std::string line;
        int n = 0;
        while (std::getline(read, line)) {
            std::istringstream readline(line);
            if (n < 8) {
                readline >> strings[n];
            } else if (n < 19) {
                readline >> doubles[n-8];
            } else {
                readline >> ints[n-19];
            }
            ++n;
        }
        read.close();
        // initialize inputs
        inputs in{
            .input_psi_file = strings[0],
            .output_psi_file = strings[1],
            .pot_file = strings[2],
            .run_type = strings[3],
            .channels = ints[5],
            .initial_pos = doubles[0],
            .final_pos = doubles[1],
            .space_grid = ints[0],
            .nx_prints = ints[3],
            .space_grid_coarse = ints[1],
            .initial_t = doubles[2],
            .final_t = doubles[3],
            .time_grid = ints[2],
            .nt_prints = ints[4],
            .potential_type = strings[4],
            .pot_read_from_file = strings[5],
            .pot_position_1 = doubles[4],
            .pot_position_2 = doubles[5],
            .strength_1 = doubles[6],
            .strength_2 = doubles[7],
            .wavepacket_type = strings[6],
            .wp_read_from_file = strings[7],
            .delta = doubles[8],
            .momentum = doubles[9],
            .wp_position = doubles[10],
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
    const std::string output = data + "/" + in.output_psi_file;
    std::ofstream wf(output, std::ios::app | std::ios::binary);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // prepare output buffer for entire set of points
    std::vector<double> buffer;
    buffer.reserve((in.time_grid / in.nt_prints + 1) * (in.space_grid / in.nx_prints) * 2);
    return {.wf = std::move(wf), .buffer = buffer};
}

// opens wavefunction output file, but with a larger buffer (by a factor of # of channels)
wfOutput openWFOutputFileCC(const inputs& in, const std::string &data) {
    const std::string output = data + "/" + in.output_psi_file;
    std::ofstream wf(output, std::ios::app | std::ios::binary);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    std::vector<double> buffer;
    buffer.reserve((in.time_grid / in.nt_prints + 1) * (in.space_grid / in.nx_prints) * 2 * in.channels);
    return {.wf = std::move(wf), .buffer = buffer};
}
