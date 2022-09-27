/**
 * Just a simple copy'n'paste template class to quickly add tests
 */

#include "Base.hh"

#include "../src/solver.hpp"

static const int TEST_ITERATIONS = 3;
using namespace std;

void test_solver() {
    const int FIELD_WIDTH = 5;
    const int FIELD_HEIGHT = 5;

    std::clock_t c_start = std::clock();

    FluidSolver *solver;

    try {
        solver = new FluidSolver(FIELD_WIDTH, FIELD_HEIGHT);
        // solver->q->u(3,3) = 2.0f;
        for ( int i = 0; i < TEST_ITERATIONS; ++i) {
            solver->fill(3);
            solver->step();
            // solver->u->print();
        }
    } catch (const std::exception &exc) {
        std::cout << "Failed to load shader: " << exc.what() << std::endl;
    }

    if(solver) delete solver;

    std::clock_t c_end = std::clock();

    long time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "cpu time: " << time_elapsed_ms << "ms" <<
              std::endl;

}

int main(int argc, const char *argv[]) {
    test_solver();
    return 0;
}