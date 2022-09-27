#include <iostream>
#include <array>
#include <future>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <unistd.h>

#include "Base.hh"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

#define COORD_TO_INDEX(ARR, X, Y, YS) ARR[X * YS + Y]

void test_vector_translation() {
    int xs = 4, ys = 4;
    std::vector<std::string> b = std::vector<std::string>(xs * ys);
    int i = 0;
    for (int y = 0; y < ys; y++) {
        for (int x = 0; x < xs; x++) {
            b[i++] = std::string("[") + std::to_string(x) + std::string(",") + std::to_string(y) + std::string("]");
        }
    }

    std::cout << "direct iteration" << std::endl;

    for (std::vector<std::string>::size_type i = 0; i != b.size(); i++) {
        std::cout << b[i] << std::endl;
    }

    std::cout << "coordinate translation" << std::endl;

    for (int y = 0; y < ys; y++) {
        for (int x = 0; x < xs; x++) {
            std::cout << "Element at " << y << ',' << x << " = " << COORD_TO_INDEX(b, x, y, ys) << std::endl;
        }
    }

}

int main(int argc, const char *argv[]) {
    test_vector_translation();
    return 0;
}
