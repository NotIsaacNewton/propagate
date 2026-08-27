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
// progress bar builder
std::string makeBar(const int progress) {
    std::string bar = "(";
    for (int i = 0; i < progress; i+=2) {
        bar += '=';
    }
    if (progress != 100) {
        bar += ">>";
    }
    for (int i = 10; i < 100 - progress; i+=2) {
        bar += ' ';
    }
    bar += ") " + std::to_string(progress) + "%";
    return bar;
}
// simple progress bar
void progressBar(const std::string& color, const int progress) {
    std::cout << color << "\r" << makeBar(progress) << RESET;
}
// reset
void reset() {
    std::cout << RESET;
}