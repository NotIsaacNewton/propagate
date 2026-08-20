//
// Created by Arian Dovald on 8/19/26.
//

#include "propagate_td.h"
#include "console_tools.h"
#include <fstream>

// creates potential operator array from data file and outputs to op
void definePotentialOperatorTD(const inputs& in, fftw_complex *op, const std::vector<double>& potential) {
    // write potential operator
    for (int i = 0; i < in.space_grid; i++) {
        const double phase = potential[i] * in.dt / 2.0;
        op[i][0] = cos(phase);
        op[i][1] = -sin(phase);
    }
}

// prepares fftw variables with RAII
fftwResources fftwPrepTD(const inputs& in, fftw_complex *psi, const std::string& data) {
    // locate fftw wisdom file
    const std::string wisdomfile = data + "/fftw_wisdom.dat";
    fftw_import_wisdom_from_filename(wisdomfile.c_str());
    // create fft and inverse fft plans with RAII
    auto fft_ptr = std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>(
            fftw_plan_dft_1d(in.space_grid, psi, psi, FFTW_FORWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    auto ifft_ptr = std::unique_ptr<std::remove_pointer_t<fftw_plan>, void(*)(fftw_plan)>(
            fftw_plan_dft_1d(in.space_grid, psi, psi, FFTW_BACKWARD, FFTW_MEASURE),
            &fftw_destroy_plan
    );
    // save fftw wisdom to file
    fftw_export_wisdom_to_filename(wisdomfile.c_str());
    // locate potential curve file
    const std::string potfile = data + "/potential.dat";
    // define only kinetic term since potential is defined at each time step
    const auto V = fftw_alloc_complex(in.space_grid);
    const auto T = fftw_alloc_complex(in.space_grid);
    defineKineticOperator(in, T, false);
    // return fftwResources
    return fftwResources{
        .fft_ptr = std::move(fft_ptr), .ifft_ptr = std::move(ifft_ptr),
        .Vp = std::unique_ptr<fftw_complex, void(*)(void*)>(V, fftw_free),
        .Tp = std::unique_ptr<fftw_complex, void(*)(void*)>(T, fftw_free)
    };
}

// prepares interpolator for potential
interpolator prepInterpolation(const inputs& in, const std::vector<std::vector<double>>& potential) {
    std::vector<double> pos_grid(in.space_grid_coarse); // stores position grid on which potential is defined
    std::vector<double> time_grid(in.time_grid_coarse); // stores time grid on which potential is defined
    const double dx = (in.final_pos-in.initial_pos)/(in.space_grid_coarse-1); // width of potential grid
    const double dt = (in.final_t-in.initial_t)/(in.time_grid_coarse-1); // width of time grid
    std::vector<spline_interp> time_interp;
    time_interp.reserve(in.space_grid_coarse);
    // write potential grids
    for (int i = 0; i < in.space_grid_coarse; i++) {
        pos_grid[i] = in.initial_pos + i*dx;
    }
    for (int i = 0; i < in.time_grid_coarse; i++) {
        time_grid[i] = in.initial_t + i*dt;
    }
    // generate interpolators for times at each position
    std::vector<double> slice(in.time_grid_coarse);
    for (int i = 0; i < in.space_grid_coarse; i++) {
        for (int j = 0; j < in.time_grid_coarse; j++) {
            slice[j] = potential[j][i];
        }
        time_interp.emplace_back(time_grid, slice);
    }
    // return interpolator struct
    return { .pos_grid = pos_grid, .time_interp = time_interp };
}

// interpolates potential across position grid at time
void interpolatePotential(const inputs& in, interpolator& interp, std::vector<double>& pot_coarse,
    std::vector<double>& potential, const int time) {
    for (int i = 0; i < in.space_grid_coarse; i++) {
        pot_coarse[i] = interp.time_interp[i].interp(time*in.dt + in.initial_t);
    }
    spline_interp pos_interp(interp.pos_grid, pot_coarse);
    for (int i = 0; i < in.space_grid; i++) {
        potential[i] = pos_interp.interp(i*in.dx + in.initial_pos);
    }
}

// propagation tick
void propTickTD(const int tick, const inputs& in, fftw_complex *psi, std::vector<double>& pot_coarse, std::vector<double>& potential,
    fftw_complex* V, const fftw_complex* T, fftw_plan fft, fftw_plan ifft, interpolator& interp, const double scale) {
    // interpolate potential
    interpolatePotential(in, interp, pot_coarse, potential, tick);
    // define potential operator
    definePotentialOperatorTD(in, V, potential);
    // apply e^(-i dt V(X))
    applyPotentialOperator(in.space_grid, psi, V);
    // execute fft plan
    fftw_execute(fft);
    // apply e^(-i dt p^2 / 2)
    applyKineticOperator(in.space_grid, psi, T);
    // execute inverse fft plan
    fftw_execute(ifft);
    // normalize fftw result (fftw uses non-normalized fft algorithm)
    scale_fftw_complex(scale, psi, in.space_grid);
    // apply e^(-i dt V(X))
    applyPotentialOperator(in.space_grid, psi, V);
}

// propagates wavefunction based on general values
// Note: FFT methods implicitly impose periodic boundary conditions
void propagateTD(const inputs& in, fftw_complex *psi, const std::string& data) {
    // prep fftw variables and plans
    auto [fft_ptr, ifft_ptr, Vp, Tp] = fftwPrepTD(in, psi, data);
    fftw_complex* V = Vp.get();
    fftw_complex* T = Tp.get();
    fftw_plan fft = fft_ptr.get();
    fftw_plan ifft = ifft_ptr.get();
    // scale for normalizing fft result
    const double scale = 1.0 / in.space_grid;
    // open output file
    const std::string output = data + "/psi_final.dat";
    std::ofstream wf(output, std::ios::app | std::ios::binary);
    if (!wf.is_open()) {
        std::cerr << "Failed to open " << output << "." << "\n";
    }
    // prepare output buffer for entire set of points
    std::vector<double> buffer;
    buffer.reserve((in.time_grid / in.nt_prints + 1) * (in.space_grid / in.nx_prints) * 2);
    // norm psi (often slightly off norm) and check norm
    std::vector<double> psi_squared(in.space_grid); // stores |psi|^2
    fftw_complex_square(psi, psi_squared); // calculates |psi|^2
    double norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared); // calculate and store norm
    scale_fftw_complex(1/sqrt(norm), psi, in.space_grid); // normalize psi
    fftw_complex_square(psi, psi_squared); // recalculate |psi|^2
    norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared); // recalculate norm
    std::cout << RED << "Check norm:\n" << RESET;
    std::print("The initial norm is {}\n", norm);
    // prepare interpolation objects
    std::vector<std::vector<double>> potential; // full coarse potential array
    readArray2D(data + "/potential.dat", potential,
        in.space_grid_coarse, in.time_grid_coarse); // reads in potential
    std::vector<double> pot_coarse(in.space_grid_coarse); // intermediate array for interpolation
    std::vector<double> pot(in.space_grid); // stores potential used to define the potential operator at each time step
    interpolator interp = prepInterpolation(in, potential); // interpolator objects
    // spacer
    spacer(RESET);
    // console output
    std::cout << "Propagation progress:\n";
    // propagation loop
    for (int t = 0; t <= in.time_grid; t++) {
        // print completion % to console
        !(t % (in.time_grid / 10)) ? std::cout << GREEN << "\r" << 100*t/in.time_grid << "%" : std::cout << RESET;
        // write lines in output file
        writeOutput(psi, t, in.space_grid, in.nx_prints, in.nt_prints, buffer);
        // propagate for one tick
        propTickTD(t, in, psi, pot_coarse, pot, V, T, fft, ifft, interp, scale);
    }
    // save buffer to output file and close the file
    wf.write(reinterpret_cast<const char*>(buffer.data()),
         static_cast<std::streamsize>(buffer.size() * sizeof(double)));
    wf.close();
    // console output
    std::cout << "\n";
    // spacer
    spacer(RESET);
    // check norm
    fftw_complex_square(psi, psi_squared);
    std::cout << RED << "Check norm:\n" << RESET;
    norm = fftw_complex_integrate(in.space_grid, in.dx, psi_squared);
    std::print("The final norm is {}\n", norm);
}

// inputs: location/of/input_file location/of/data_directory
int main(const int argc, const char* argv[]) {
    if (argc != 3) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}[location/of/input_file] [location/of/data_directory]\n", GREEN);
        spacerFancy(RED);
        return 1;
    }

    // record the start time
    auto start = std::chrono::steady_clock::now();

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];
    const std::string psifile = data + "/psi_initial.dat";
    const std::string psiout = data + "/psi_final.dat";

    // spacer
    spacerChunky(BLUE);
    std::print("\n");

    // read input file
    const inputs in = readInputs(inputfile);

    // warning if odd number of gridpoints
    if (in.space_grid % 2 != 0) {
        spacerFancy(YELLOW);
        std::cerr << YELLOW << "WARNING: Odd number of gridpoints may produce asymmetric momentum grid.\n" << RESET;
        spacerFancy(YELLOW);
    }
    // warning if large dt
    if (in.time_grid < 1000) {
        spacerFancy(YELLOW);
        std::cerr << YELLOW << "WARNING: Large dt may fail to resolve rapid dynamics. Use dt < pi/E_max.\n" << RESET;
        spacerFancy(YELLOW);
    }

    // clear output file
    std::ofstream clearout;
    clearout.open(psiout);
    clearout.close();

    // allocate wavefunction with RAII
    auto psi = fftw_alloc_complex(in.space_grid);
    std::unique_ptr<fftw_complex, void(*)(void*)> psip{psi, fftw_free};
    // get wavefunction from psifile and save it to psi
    fftw_complex_array_from_file(psifile, psi, in.space_grid);

    // spacer
    spacerThick(RESET);

    // propagate wf
    propagateTD(in, psi, data);

    // spacer
    spacerThick(RESET);

    // record end time and duration
    auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    std::print("Done :)\n\n");
    std::print("{}Execution time: {} seconds\n", BLUE, sec.count());

    // spacer
    spacerChunky(BLUE);

    return 0;
}
