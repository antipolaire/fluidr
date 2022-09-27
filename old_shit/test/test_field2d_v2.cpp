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

#include "../src/field2d_v2.hpp"

using namespace Eigen;
using namespace std;

void test_field2d() {
    Field2DV2<float> *u = new Field2DV2<float>(10, 10);

    cout << "Test print Field2D<Vector2f>" << std::endl;
    cout << std::endl;

    cout << u << std::endl;

    cout << std::endl;
    cout << "Test print Field2D<Vector2f> by lines" << std::endl;
    cout << std::endl;

    u->forEach([](auto &_u, auto &_v, int x, int y) {
        cout << "F[" << x << "][" << y << "] = (" << _u << "," << _v << ")" << endl;
    });

    cout << std::endl;
    cout << "Test init Field2D<Vector2f> with coordinates" << std::endl;
    cout << std::endl;

    u->forEach([](auto &u_xy, auto &v_xy, int x, int y) {
        u_xy = x;
        v_xy = y;
    });

    cout << std::endl;
    cout << "Test print Field2D<Vector2f> with coordinates" << std::endl;
    cout << std::endl;


    u->forEach([](auto &_u, auto &_v, int x, int y) {
        cout << "F[" << x << "][" << y << "] = (" << _u << "," << _v << ")" << endl;
    });

    cout << std::endl;
    cout << "Test init Field2D<float>" << std::endl;
    cout << std::endl;

    Field2DV2<float> *u2 = new Field2DV2<float>(10, 10);

    // u2.print(cout);

    cout << std::endl;
    cout << "Test print Field2D<float> by lines" << std::endl;
    cout << std::endl;

    u->forEach([](auto &_u, auto &_v, int x, int y) {
        cout << "F[" << x << "][" << y << "] = (" << _u << ", " <<  _v << ")" << endl;
    });


    cout << std::endl;
    cout << "Test init Field2D<float> coordinate of x" << std::endl;
    cout << std::endl;

    u->forEach([](auto &u_xy, auto &v_xy, int x, int y) {
        u_xy = x;
    });

    cout << std::endl;
    cout << "Test print Field2D<float> with new values" << std::endl;
    cout << std::endl;

    u->forEach([](auto &_u, auto &_v, int x, int y) {
        cout << "F[" << x << "][" << y << "] = ("  << _u << ", " <<  _v <<  ")" << endl;
    });

    delete u;
    delete u2;
}

void test_setter_vs_getter(){
    Field2DV2<float> *field = new Field2DV2<float>(10, 10);

    for ( int i = 0; i < 10*10; ++i ) {
        field->u(i) = (float)i;
        field->v(i) = (float)i;
    }

    for ( int i = 0; i < 10*10; ++i ) {
        ALEPH_ASSERT_EQUAL(field->u(i), (float)i)
        ALEPH_ASSERT_EQUAL(field->v(i), (float)i)
    }

    field->print();

    delete field;
}

int main(int argc, const char *argv[]) {
    test_field2d();
    test_setter_vs_getter();
    return 0;
}
