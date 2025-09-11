//
// Created by Arian Dovald.
//

#include "filetools.h"
#include "propagate.h"
#include "potentials.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"

// TODO: more safety checks and error paths (try <expected>)

// inputs: location/of/input_file location/of/data_directory
int main(int argc, const char* argv[]) {
    if (argc != 3) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}location/of/input_file location/of/data_directory\n", GREEN);
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

    // clear output file
    std::ofstream clearout;
    clearout.open(psiout);
    clearout.close();

    // allocate wavefunction with RAII
    auto psi = fftw_alloc_complex(in.space_grid);
    std::unique_ptr<fftw_complex, void(*)(void*)> psip{psi, fftw_free};
    // get wavefunction from psifile
    fftw_complex_array_from_file(psifile, psi);

    // spacer
    spacerThick(RESET);

    // propagate wf
    propagate(in, psi, data);

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
