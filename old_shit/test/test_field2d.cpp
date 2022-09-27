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

#include "../src/field2d.hpp"

using namespace Eigen;
using namespace std;

void test_field2d() {
    Field2D<Vector2f> *u = new Field2D<Vector2f>(10, 10);

    cout << "Test print Field2D<Vector2f>" << std::endl;
    cout << std::endl;

    cout << u << std::endl;

    cout << std::endl;
    cout << "Test print Field2D<Vector2f> by lines" << std::endl;
    cout << std::endl;

    u->forEach([](auto &cell, int x, int y) {
        cout << "F[" << x << "][" << y << "] = (" << cell.x() << "," << cell.y() << ")" << endl;
    });

    cout << std::endl;
    cout << "Test init Field2D<Vector2f> with coordinates" << std::endl;
    cout << std::endl;

    u->forEach([](auto &cell, int x, int y) {
        cell << x, y;
    });

    cout << std::endl;
    cout << "Test print Field2D<Vector2f> with coordinates" << std::endl;
    cout << std::endl;


    u->forEach([](auto &cell, int x, int y) {
        cout << "F[" << x << "][" << y << "] = (" << cell.x() << "," << cell.y() << ")" << endl;
    });

    cout << std::endl;
    cout << "Test init Field2D<float>" << std::endl;
    cout << std::endl;

    Field2D<float> *u2 = new Field2D<float>(10, 10);

    // u2.print(cout);

    cout << std::endl;
    cout << "Test print Field2D<float> by lines" << std::endl;
    cout << std::endl;

    u2->forEach([](auto &cell, int x, int y) {
        cout << "F[" << x << "][" << y << "] = (" << cell << ")" << endl;
    });


    cout << std::endl;
    cout << "Test init Field2D<float> coordinate of x" << std::endl;
    cout << std::endl;

    u2->forEach([](auto &cell, int x, int y) {
        cell = x;
    });

    cout << std::endl;
    cout << "Test print Field2D<float> with new values" << std::endl;
    cout << std::endl;

    u2->forEach([](auto &cell, int x, int y) {
        cout << "F[" << x << "][" << y << "] = (" << cell << ")" << endl;
    });

    delete u;
    delete u2;
}

int main(int argc, const char *argv[]) {
    test_field2d();
    return 0;
}
