#include <iostream>
#include <array>
#include <future>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <unistd.h>

#include "Base.hh"

#include "../src/solver.hpp"

using namespace Eigen;
using namespace std;

void test_modular_solver() {

    auto *solver = new FluidSolver(50, 50);
    // solver->u->print();

    std::clock_t c_start = std::clock();

    for (int i = 0; i < 1; ++i) {
        solver->step();
        //solver->u->print();
    }
    std::clock_t c_end = std::clock();

    unsigned long time_elapsed_ms = 1000 * (c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "cpu time: ";
    std::cout << time_elapsed_ms;
    std::cout << "ms";
    std::cout << std::endl;

    delete solver;
}

int main(int argc, const char *argv[]) {
    test_modular_solver();
    return 0;
}
