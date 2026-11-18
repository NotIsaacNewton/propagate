//
// Created by Arian Dovald on 11/18/26.
//

#include "propagate_td_cc_lr_yoshi.h"
#include "fftw_complex_tools.h"
#include "console_tools.h"
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <print>
#include <algorithm>
#include <variant>
#include <omp.h>
#include "../potentials/potentials_td_cc.h"
#include "../potentials/potentials_td.h"


