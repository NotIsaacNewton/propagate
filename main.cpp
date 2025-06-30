//#include <cstdlib>
#include "filetools.h"
#include "propagate.h"
#include "fftw_complex_tools.h"
#include "wavepackets_and_potentials.h"

// ANSI escape codes for colors
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
//#define BLUE    "\033[34m"

// TODO: create header that automatically selects between ArmPL and MKL (urgent, to run on beocat)
// TODO: make everything more general and flexible for use in future projects

int main() {
    // record the start time
    auto start = std::chrono::high_resolution_clock::now();

    // file locations
    std::string inputfile = "/Users/ariandovald/CLionProjects/propagate/inputs";
    std::string psifile = "/Users/ariandovald/CLionProjects/propagate/psi_initial.dat";
    std::string psisquared = "/Users/ariandovald/CLionProjects/propagate/psi_squared_initial.dat";
    std::string psiout = "/Users/ariandovald/CLionProjects/propagate/psi_final.dat";
    std::string psisquaredout = "/Users/ariandovald/CLionProjects/propagate/psi_squared_final.dat";

    // spacer
    std::cout << YELLOW << "--------------------------------------------------------------------------------\n" << RESET;

    // space stuff
    double x0;
    double xf;
    int ndx;
    // time stuff
    double t0;
    double tf;
    int ndt;

    // read input file
    readInputs(x0, xf, ndx, t0, tf, ndt, inputfile);
    // space-time grid-widths
    const double dx = (xf - x0)/ndx;
    const double dt = (tf - t0)/ndt;

    // clear output file
    std::ofstream potwrite;
    potwrite.open(psiout);
    potwrite.close();

    // spacer
    std::cout << YELLOW << "--------------------------------------------------------------------------------\n" << RESET;

    // write function to array, write array to file
    fftw_complex psi[ndx];
    fftw_complex_func_to_array(x0,xf,dx,gaussianWP,psi);
    writeComplex1D(x0, xf, dx, psifile, psi);

    // save initial wf-squared to file
    double psi_squared[ndx];
    fftw_complex_square(ndx, psi, psi_squared);
    writeArray1D(x0, xf, dx, psisquared, psi_squared);

    // spacer
    std::cout << YELLOW << "--------------------------------------------------------------------------------\n" << RESET;

    // check norm
    std::cout << RED "Check norm:\n" << RESET;
    double norm = array_integrate(ndx, dx, psi_squared);
    std::print("The initial norm is {}\n", norm);

    // spacer
    std::cout << YELLOW << "--------------------------------------------------------------------------------\n" << RESET;

    // propagate wf
    std::cout << "Propogation progress:\n";
    propagate(ndx, dx, x0, psi, step, dt, ndt, psiout);

    // spacer
    std::cout << YELLOW << "--------------------------------------------------------------------------------\n" << RESET;

    // save final wf-squared to file
    fftw_complex_square(ndx, psi, psi_squared);
    writeArray1D(x0, xf, dx, psisquaredout, psi_squared);

    // spacer
    std::cout << YELLOW << "--------------------------------------------------------------------------------\n" << RESET;

    // check norm
    std::cout << RED "Check norm:\n" << RESET;
    norm = array_integrate(ndx, dx, psi_squared);
    std::print("The final norm is {}\n", norm);

    // spacer
    std::cout << YELLOW << "--------------------------------------------------------------------------------\n" << RESET;

    // record end time and duration
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time = double(duration.count())/1000000;
    std::cout << "Done :) Execution time: " << time << " seconds\n";

    // TODO: would be nice to pipe commands to mathematica or gnuplot

    return 0;
}
