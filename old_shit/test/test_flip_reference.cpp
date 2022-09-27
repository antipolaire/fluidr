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


void assert_equals_vector_value(Vector2d *&v, const Vector2d &vv) {
    ALEPH_ASSERT_EQUAL(v->x(), vv.x())
    ALEPH_ASSERT_EQUAL(v->y(), vv.y())
}

/**
 * Idea:
 * 1. Have data-structures 'a' and 'b'.
 * 2. Copy content of 'a' to 'b'.
 * 3. Modify 'b'
 * 4. Flip address of 'a' and 'b'
 */
void test_flip_reference() {
    Vector2d *a = new Vector2d();
    Vector2d *b = new Vector2d();

    (*a) << 1.0, 1.0;
    (*b) << 2.0, 2.0;

    assert_equals_vector_value(a, Vector2d(1.0, 1.0));
    assert_equals_vector_value(b, Vector2d(2.0, 2.0));

    memcpy(b->data(), a->data(), a->size() * sizeof(double));

    (*b) << 3.0, 3.0;

    auto tmp = a;
    a = b;
    b = tmp;


    assert_equals_vector_value(a, Vector2d(3.0, 3.0));
    assert_equals_vector_value(b, Vector2d(1.0, 1.0));

    delete a;
    delete b;
}

void test_flip_reference_field() {
    int dimension_x = 3, dimension_y = 3;

    Field2D<Vector2d> *a = new Field2D<Vector2d>(dimension_x, dimension_y);

    Field2D<Vector2d> *b = new Field2D<Vector2d>(dimension_x, dimension_y);

    a->forEach([](auto &cell, int x, int y) {
        cell << 0.0, 0.0;
    });


    b->forEach([](auto &cell, int x, int y) {
        cell << 10.0, 10.0;
    });

    b->copyTo(a);

    a->forEach([&](auto &cell, int x, int y) {
        auto a_value = (*a)(x, y);
        auto b_value = (*b)(x, y);
        ALEPH_ASSERT_EQUAL(a_value.x(), b_value.x())
        ALEPH_ASSERT_EQUAL(a_value.y(), b_value.y())
    });

    delete a;
    delete b;
}

int main(int argc, const char *argv[]) {
    test_flip_reference();
    test_flip_reference_field();
    return 0;
}
