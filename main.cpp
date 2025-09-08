//
// Created by Arian Dovald.
//

#include "filetools.h"
#include "propagate.h"
#include "potentials.h"
#include "fftw_complex_tools.h"
#include "wavepackets.h"
#include "console_tools.h"

// TODO:
//  separate the propagation code from potential creation and wavefunction creation
//  create header that automatically selects between ArmPL and MKL
//  <expected> for error paths
int main(int argc, const char* argv[]) {
    // record the start time
    auto start = std::chrono::steady_clock::now();

    // file locations
    std::string inputfile = argv[1];
    std::string data = argv[2];
    std::string psifile = data+"/psi_initial.dat";
    std::string psiout = data+"/psi_final.dat";

    // spacer
    spacerChunky(BLUE);
    std::print("\n");

    // read input file
    inputs in = readInputs(inputfile);

    // clear output file
    std::ofstream potwrite;
    potwrite.open(psiout);
    potwrite.close();

    // allocate wavefunction with RAII
    auto psi = fftw_alloc_complex(in.space_grid);
    std::unique_ptr<fftw_complex, void(*)(void*)> psip{psi, fftw_free};
    // write wavefunction to array, save array to file
    fftw_complex_func_to_array(in.initial_pos,in.final_pos,in.dx,
        gaussianWP(0.5, 10, 3), psi);
    fftw_complex_array_to_file(in.initial_pos, in.final_pos, in.dx,
        psifile, psi);
    // NOTE: replace above calls with fftw_complex_array_from_file to read wf from an input file

    // spacer
    spacerThick(RESET);

    // propagate wf
    propagate(in, psi, sho(5), data);

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
