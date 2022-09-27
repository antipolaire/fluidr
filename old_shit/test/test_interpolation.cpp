/**
 * Just a simple copy'n'paste template class to quickly add tests
 */

#include "Base.hh"
#include "../src/math.hpp"
#include "../src/field2d.hpp"

using namespace std;
//
//void assert1DInterpolation2ndOrder(const double *u, const function<double(double)> &kernel, double offset,
//                                   double expectation);
//
//void assert1DInterpolation4thOrder(double *u, const function<double(double)> &kernel, double offset,
//                                   double expectation, const double radius);
//
//void assert2DInterpolation4thOrder(const function<double(double)> &kernel, double kernel_radius, Field2D<Vector2d> *&u,
//                                   Field2D<Vector2d> *&q, Field2D<Vector2d> *&expected, double alpha,
//                                   double beta);
//
//using namespace std;
//
//void assert1DInterpolation2ndOrder(const double *u, const function<double(double)> &kernel, double offset,
//                                   double expected) {
//    double actual = fluid::util::math::interpolation::general::apply_2nd_order(u[0], u[1], kernel, offset);
//    std::cout << "offset=" << offset << std::endl
//              << "expected=" << expected << std::endl
//              << "actual=" << actual << std::endl << std::endl;
//    ALEPH_ASSERT_EQUAL(actual, expected)
//}
//
//void assert1DInterpolation4thOrder(double *u, const function<double(double)> &kernel, double offset,
//                                   double expected, const double radius = 1.0) {
//    double actual = fluid::util::math::interpolation::general::apply_4th_order(u[0], u[1], u[2], u[3],
//                                                                               kernel, offset, radius);
//    // std::cout << offset << " => expected=" << expected << " : actual:" << actual << std::endl;
//    ALEPH_ASSERT_EQUAL(actual, expected)
//}
//
//
//void assert2DInterpolation4thOrder(const function<double(double)> &kernel, double kernel_radius, Field2D<Vector2d> *&u,
//                                   Field2D<Vector2d> *&q, Field2D<Vector2d> *&expected, double alpha,
//                                   double beta) {
//
//    function<double(int, int)> f_x = [&q](int x, int y) {
//        return (*q)(x, y).x();
//    };
//
//    function<double(int, int)> f_y = [&q](int x, int y) {
//        return (*q)(x, y).y();
//    };
//
//    cout << "actual[before]: " << u << endl;
//
//    u->forEach([&](auto &cell, int x, int y) {
//        cell.x() = fluid::util::math::interpolation::general::bi_apply_4th_order(-1, 0, 1, 2, x, y, f_x, kernel,
//                                                                                 alpha,
//                                                                                 beta,
//                                                                                 kernel_radius);
//
//        cell.y() = fluid::util::math::interpolation::general::bi_apply_4th_order(-1, 0, 1, 2, x, y, f_y, kernel,
//                                                                                 alpha,
//                                                                                 beta,
//                                                                                 kernel_radius);
//    });
//
//    cout << "actual[after]: " << u << endl;
//
//    // cout << "expected: " << expected << endl;
//
//    for (int i = 0; i < 5; i++) {
//        for (int j = 0; j < 5; j++) {
////            ALEPH_ASSERT_EQUAL(v_expected(i, j).x(), v_actual(i, j).x())
////            ALEPH_ASSERT_EQUAL(v_expected(i, j).y(), v_actual(i, j).y())
//        }
//    }
//}
//
//
//void test_interp_1d_2nd_order_1() {
//
//    double u[2];
//
//    std::function<double(double)> kernel = fluid::util::math::interpolation::general::kernel_bilinear;
//
//    std::cout << std::endl << "test 1d 2nd order interpolation - range [0.0,1.0]:" << std::endl << std::endl;
//
//    u[0] = 0.0;
//    u[1] = 1.0;
//
//    assert1DInterpolation2ndOrder(u, kernel, -1.0, 0.0);
//    assert1DInterpolation2ndOrder(u, kernel, -0.5, 0.0);
//    assert1DInterpolation2ndOrder(u, kernel, 0.0, 0.0);
//    assert1DInterpolation2ndOrder(u, kernel, 0.25, 0.25);
//    assert1DInterpolation2ndOrder(u, kernel, 0.5, 0.5);
//    assert1DInterpolation2ndOrder(u, kernel, 0.75, 0.75);
//    assert1DInterpolation2ndOrder(u, kernel, 1.0, 1.0);
//    assert1DInterpolation2ndOrder(u, kernel, 1.5, 1.0);
//    assert1DInterpolation2ndOrder(u, kernel, 2.0, 1.0);
//
//}
//
//void test_interp_1d_2nd_order_2() {
//
//    double u[2];
//
//    std::function<double(double)> kernel = fluid::util::math::interpolation::general::kernel_bilinear;
//
//    std::cout << std::endl << "test 1d 2nd order interpolation - range [1.0,2.0]:" << std::endl << std::endl;
//
//
//    u[0] = 1.0;
//    u[1] = 2.0;
//
//    assert1DInterpolation2ndOrder(u, kernel, -1.0, 0.0);
//    assert1DInterpolation2ndOrder(u, kernel, -0.5, 0.5);
//    assert1DInterpolation2ndOrder(u, kernel, 0.0, 1.0);
//    assert1DInterpolation2ndOrder(u, kernel, 0.25, 1.25);
//    assert1DInterpolation2ndOrder(u, kernel, 0.5, 1.5);
//    assert1DInterpolation2ndOrder(u, kernel, 0.75, 1.75);
//    assert1DInterpolation2ndOrder(u, kernel, 1.0, 2.0);
//    assert1DInterpolation2ndOrder(u, kernel, 1.5, 2.0);
//    assert1DInterpolation2ndOrder(u, kernel, 2.0, 2.0);
//}
//
//void test_interp_1d_2nd_order_3() {
//
//    double u[2];
//
//    std::function<double(double)> kernel = fluid::util::math::interpolation::general::kernel_bilinear;
//
//    std::cout << std::endl << "test 1d 2nd order interpolation - range [-1.0,0.0]:" << std::endl << std::endl;
//
//    u[0] = -1.0;
//    u[1] = 0.0;
//
//    assert1DInterpolation2ndOrder(u, kernel, -1.0, -2.0);
//    assert1DInterpolation2ndOrder(u, kernel, -0.5, -1.5);
//    assert1DInterpolation2ndOrder(u, kernel, 0.0, -1.0);
//    assert1DInterpolation2ndOrder(u, kernel, 0.25, -0.75);
//    assert1DInterpolation2ndOrder(u, kernel, 0.5, -0.5);
//    assert1DInterpolation2ndOrder(u, kernel, 0.75, -0.25);
//    assert1DInterpolation2ndOrder(u, kernel, 1.0, 0.0);
//    assert1DInterpolation2ndOrder(u, kernel, 1.5, 0.0);
//    assert1DInterpolation2ndOrder(u, kernel, 2.0, 0.0);
//
//}
//
//void test_interp_1d_4th_order() {
//
//    std::function<double(double)> kernel = fluid::util::math::interpolation::general::kernel_bilinear;
//    const double radius = 2.0;
//
//    double u[4];
//
//    u[0] = 0.0;
//    u[1] = 1.0;
//    u[2] = 2.0;
//    u[3] = 3.0;
//
//    assert1DInterpolation4thOrder(u, kernel, -2.0, 0.0, radius);
//    assert1DInterpolation4thOrder(u, kernel, -1.5, 0.0, radius);
//    assert1DInterpolation4thOrder(u, kernel, -1.0, 0.0, radius);
//    assert1DInterpolation4thOrder(u, kernel, -0.5, 0.5, radius);
//    assert1DInterpolation4thOrder(u, kernel, 0.0, 1.0, radius);
//    assert1DInterpolation4thOrder(u, kernel, 0.25, 1.25, radius);
//    assert1DInterpolation4thOrder(u, kernel, 0.5, 1.5, radius);
//    assert1DInterpolation4thOrder(u, kernel, 0.75, 1.75, radius);
//    assert1DInterpolation4thOrder(u, kernel, 1.0, 2.0, radius);
//    assert1DInterpolation4thOrder(u, kernel, 1.5, 2.5, radius);
//    assert1DInterpolation4thOrder(u, kernel, 2.0, 3.0, radius);
//    assert1DInterpolation4thOrder(u, kernel, 2.5, 3.0, radius);
//    assert1DInterpolation4thOrder(u, kernel, 3.0, 3.0, radius);
//
//}
//
//void test_bridson_interp() {
//
//    std::function<double(double)> kernel = fluid::util::math::interpolation::general::kernel_bilinear;
//
//    double u[4][4];
//
//    u[0][0] = 0.0;
//    u[0][1] = 1.0;
//    u[0][2] = 2.0;
//    u[0][3] = 3.0;
//
//    u[1][0] = 1.0;
//    u[1][1] = 1.0;
//    u[1][2] = 2.0;
//    u[1][3] = 3.0;
//
//    u[2][0] = 2.0;
//    u[2][1] = 2.0;
//    u[2][2] = 2.0;
//    u[2][3] = 3.0;
//
//    u[3][0] = 3.0;
//    u[3][1] = 3.0;
//    u[3][2] = 3.0;
//    u[3][3] = 3.0;
//
//    std::function<double(int, int)> f = [&u](int x, int y) {
//        auto _x = std::clamp(x, 0, 3);
//        auto _y = std::clamp(y, 0, 3);
//        return u[_x][_y];
//    };
//
//    double alpha = 0.0;
//    double beta = 0.0;
//
//    for (int x = 0; x < 4; x++) {
//        for (int y = 0; y < 4; y++) {
//            auto value = fluid::util::math::interpolation::general::bi_apply_4th_order(-1.0, 0.0, 1.0, 2.0, x, y, f,
//                                                                                       kernel, alpha, beta);
//            // double value = fluid::util::math::interpolation::bridson(f, x, y, 3.0, 3.0);
//            std::cout << value << " ";
//        }
//        std::cout << std::endl;
//    }
//
//}
//
//void test_interp_2d_4th_order() {
//
//
//    std::function<double(double)> kernel = fluid::util::math::interpolation::general::kernel_bilinear;
//
//    auto *u = new Field2D<Vector2d>(5, 5);
//    auto *q = new Field2D<Vector2d>(5, 5);
//    auto *expected = new Field2D<Vector2d>(5, 5);
//
//    (*q)(2, 2) = Vector2d(1.0, 1.0);
///*
//    (*expected)(2, 2) = Vector2d(0.5, 0.5);
//    (*expected)(2, 3) = Vector2d(0.0, 0.0);
//    (*expected)(2, 1) = Vector2d(0.0, 0.0);
//    (*expected)(3, 2) = Vector2d(0.5, 0.5);
//    (*expected)(1, 2) = Vector2d(0.0, 0.0);
//*/
//
//    for (int i = 0; i < 1; i++) {
//        assert2DInterpolation4thOrder(kernel, 1.0, u, q, expected, 1.0, 0.9);
//        u->copyTo(q);
//    }
//
//    delete u;
//    delete q;
//    delete expected;
//}

int main(int argc, const char *argv[]) {


    std::cout << fluid::util::math::interpolation::lerp(0.0, 1.0, 1.0) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(0.0, 1.0, 0.5) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(0.0, 1.0, 0.0) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(0.0, 1.0, -.5) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(0.0, 1.0, -1.0) << std::endl;

    std::cout << fluid::util::math::interpolation::lerp(-1.0, 0.0, 1.0) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(-1.0, 0.0, 0.5) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(-1.0, 0.0, 0.0) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(-1.0, 0.0, -.5) << std::endl;
    std::cout << fluid::util::math::interpolation::lerp(-1.0, 0.0, -1.0) << std::endl;


    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 1.0, 0.0) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 0.5, 0.0) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 0.0, 0.0) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, -0.5, 0.0) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, -1.0, 0.0) << std::endl;

    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 0.0, 1.0) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 0.0, 0.5) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 0.0, 0.0) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 0.0, -0.5) << std::endl;
    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 1.0, 1.0, 0.0, -1.0) << std::endl;

    std::cout << fluid::util::math::interpolation::bilerp(0.0, 0.0, 2.0, 1.0, 0.5, 0.5) << std::endl;

    // test_interp_1d_2nd_order_1();
    //test_interp_1d_2nd_order_2();
    //test_interp_1d_2nd_order_3();
    //test_interp_1d_4th_order();
    //test_bridson_interp();
    //test_interp_2d_4th_order();
    return 0;
}