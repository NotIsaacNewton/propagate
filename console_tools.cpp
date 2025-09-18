//
// Created by Arian Dovald on 9/18/25.
//

#include <iostream>
#include "console_tools.h"

// spacer
void spacer(const std::string& color) {
    std::cout << color << "--------------------------------------------------------------------------------\n" << RESET;
}
// thick spacer
void spacerThick(const std::string& color) {
    std::cout << color << "================================================================================\n" << RESET;
}
// chunky spacer
void spacerChunky(const std::string& color) {
    std::cout << color << "################################################################################\n" << RESET;
}
// fancy spacer
void spacerFancy(const std::string& color) {
    std::cout << color << "<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>\n" << RESET;
}
// chain spacer
void spacerChain(const std::string& color) {
    std::cout << color << "<=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=><=>\n" << RESET;
}
