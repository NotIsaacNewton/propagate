//
// Created by Arian Dovald on 7/9/26.
//

#include "interpolate_1d.h"
#include "console_tools.h"
#include <iostream>
#include "filetools.h"

// inputs: location/of/input_output_file gridpoints_in gridpoints_out initial_abscissa final_abscissa
int interpolate(const int argc, char *argv[]) {
    if(argc < 6) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}location/of/input_file gridpoints_in gridpoints_out initial_abscissa final_abscissa\n",
            GREEN);
        spacerFancy(RED);
        return 1;
    }

    // spacer
    spacerChunky(BLUE);
    std::print("\n");

    // file locations
    const std::string datafile = argv[1];
    const int gridpoints_in = std::stoi(argv[2]);
    const int gridpoints = std::stoi(argv[3]);
    const double initial_x = std::stod(argv[4]);
    const double final_x = std::stod(argv[5]);

    // where the data goes
    std::vector<double> data_in(gridpoints_in);
    std::vector<double> data(gridpoints);

    // read input
    std::print("Reading {}\n", datafile);
    readArray1D(datafile, data_in);

    // write initial grid
    std::vector<double> grid(gridpoints_in);
    const double dx_in = (final_x - initial_x)/(gridpoints_in - 1); // width of initial grid
    for (int i = 0; i < gridpoints_in; i++) {
        grid[i] = initial_x + i*dx_in;
    }

    // interpolate
    spline_interp interpolator(grid, data_in); // spline interpolation object
    const double dx = (final_x - initial_x)/(gridpoints - 1);
    // write interpolated data
    for (int i = 0; i < gridpoints; i++) {
        data[i] = interpolator.interp(initial_x + i*dx);
    }

    // send interpolated data to file
    writeArray1D(initial_x, final_x, dx, gridpoints, datafile, data);

    // spacer
    spacer(RESET);

    // console output
    std::print("{}Interpolation completed!\n\n", GREEN);
    spacerChunky(BLUE);

    return 0;
}
