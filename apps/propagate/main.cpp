//
// Created by Arian Dovald on 9/2/26.
//

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <print>
#include "file_tools.h"
#include "console_tools.h"
#include "propagate.h"
#include "propagate_td.h"
#include "propagate_td_cc_seq.h"
#include "propagate_td_cc_lr.h"
#include "propagate_td_cc_dense.h"
#include "propagate_td_cc_lr_bm02.h"

// inputs: location/of/input_file location/of/data_directory
int main(const int argc, const char* argv[]) {
    if (argc != 3) {
        spacerFancy(RED);
        std::cerr << RED << "Error: improper inputs.\n";
        std::print("{}[location/of/input_file] [location/of/data_directory]\n", GREEN);
        spacerFancy(RED);
        return 1;
    }

    // file locations
    const std::string inputfile = argv[1];
    const std::string data = argv[2];

    // introduction
    std::print("\n{}propagate\n", BLUE);

    // record the start time
    auto start = std::chrono::steady_clock::now();

    // spacer
    spacerChunky(BLUE);

    // read input file
    const inputs in(inputfile);

    // spacer
    spacer(RESET);

    // confirm run type in console
    std::print("{}——> using {}{}{} propagation{}\n", GRAY, MAGENTA, in.run_type, GRAY, RESET);

    // psi output file location
    const std::string psiout = data + "/" + in.output_psi_file;

    // imaginary propagation?
    const bool imProp = in.run_type == "IM";
    // TD propagation?
    const bool propTD = in.run_type == "TD";
    // TDCC real-space two propagation?
    const bool propTDCCSeq = in.run_type == "TDCCSeq";
    // TDCC dense-diagonalization propagation?
    const bool propTDCCDense = in.run_type == "TDCCDense";
    // TDCC low-rank propagation?
    const bool propTDCCLR = in.run_type == "TDCCLR";
    // TDCC low-rank propagation split by BM02?
    const bool propTDCCLRBM02 = in.run_type == "TDCCLRBM02";

    // warning if odd number of gridpoints
    if (in.space_grid % 2 != 0) {
        spacerFancy(RED);
        std::cerr << RED << "ERROR: Odd number of gridpoints produces asymmetric momentum grid.\n" << RESET;
        spacerFancy(RED);
        return 1;
    }
    // warning if large dt
    if (in.time_grid < 1000) {
        spacerFancy(YELLOW);
        std::cerr << YELLOW << "WARNING: Large dt may fail to resolve rapid dynamics. Try dt < pi/E_max.\n" << RESET;
        spacerFancy(YELLOW);
    }

    // clear output file
    std::ofstream clearout;
    clearout.open(psiout);
    clearout.close();

    // spacer
    spacerThick(RESET);

    // propagate wf
    if (propTD) {
        propagateTD(in, data);
    } else if (propTDCCSeq) {
        propagateTDCCSeq(in, data);
    } else if (propTDCCDense) {
        propagateTDCCDense(in, data);
    } else if (propTDCCLR) {
        propagateTDCCLR(in, data);
    } else if (propTDCCLRBM02) {
        propagateTDCCLRBM02(in, data);
    } else {
        propagate(in, data, imProp);
    }

    // spacer
    spacerThick(RESET);

    // record end time and duration
    auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> sec = end - start;
    std::print("Done :)\n");
    spacerChunky(BLUE);
    std::print("{}Execution time: {} seconds\n\n", BLUE, sec.count());

    return 0;
}
