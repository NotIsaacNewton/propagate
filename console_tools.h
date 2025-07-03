//
// Created by Arian Dovald on 7/1/25.
//

#ifndef CONSOLE_TOOLS_H
#define CONSOLE_TOOLS_H

#include <iostream>
#include <string>

// ANSI escape codes for colors
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"

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
#endif //CONSOLE_TOOLS_H
