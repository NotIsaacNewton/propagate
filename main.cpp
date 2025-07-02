//#include <cstdlib>
#include "filetools.h"
#include "propagate.h"
#include "fftw_complex_tools.h"
#include "wavepackets_and_potentials.h"
#include "console_tools.h"

// TODO: create header that automatically selects between ArmPL and MKL (urgent, to run on beocat)
// TODO: make everything more general and flexible for use in future projects

int main() {
    // record the start time
    auto start = std::chrono::high_resolution_clock::now();

    // file locations
    std::string inputfile = "/Users/ariandovald/CLionProjects/propagate/inputs";
    std::string psifile = "/Users/ariandovald/CLionProjects/propagate/data/psi_initial.dat";
    std::string psisquared = "/Users/ariandovald/CLionProjects/propagate/data/psi_squared_initial.dat";
    std::string psiout = "/Users/ariandovald/CLionProjects/propagate/data/psi_final.dat";

    // spacer
    spacerChunky(BLUE);
    std::cout << std::endl;

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

    // TODO: make it possible to choose initial wf source and parameters from input file
    // write wavefunction to array, save array to file
    fftw_complex psi[ndx];
    fftw_complex_func_to_array(x0,xf,dx,gaussianWP,psi);
    fftw_complex_array_to_file(x0, xf, dx, psifile, psi);
    // NOTE: replace above calls with fftw_complex_array_from_file to read wf from an input file

    // spacer
    spacerThick(RESET);

    // TODO: make it possible to choose potential source and parameters from input file
    // propagate wf
    propagate(ndx, dx, x0, psi, barrier(0, 0.2, 200),
              dt, ndt, psiout);

    // spacer
    spacer(RESET);

    // record end time and duration
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time = double(duration.count())/1000000;
    std::cout << "Done :)\n\n";
    std::cout << BLUE << "Execution time: " << time << " seconds\n";

    // spacer
    spacerChunky(BLUE);

    // TODO: would be nice to pipe commands to mathematica or gnuplot

    return 0;
}
