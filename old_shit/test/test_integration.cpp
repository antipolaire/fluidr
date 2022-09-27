/**
 * Just a simple copy'n'paste template class to quickly add tests
 */

#include "Base.hh"

#include "../src/bnas.h"
#include "../src/math.hpp"

using namespace std;

double _rk4(const std::function<double(double, double)> &f, double x, double y, double h) {
    auto k1 = h * f(x, y);
    auto k2 = h * f(x + h / 2.0, y + k1 / 2.0);
    auto k3 = h * f(x + h / 2.0, y + k2 / 2.0);
    auto k4 = h * f(x + h, y + k3);
    return y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

void test_integrate_point() {

    // Idea: approximating f only knowing f'

    // f
    std::function<double(double)> f = [](double x) {
        return (x * x + 4) * (x * x + 4) / 16;
    };

    // f'
    std::function<double(double, double)> f_ = [&](double x, double y) {
        return x * std::sqrt(y);
    };

    double dt = 0.5;

    double y = 1.0;
    for (double x = 0.0; x < 10.0; x += dt) {
        y = _rk4(f_, x, y, dt);
        std::cout << "f  (" << x + dt << ") = " << f(x + dt) << std::endl
                  << "~f (" << x + dt << ") = " << y << std::endl
                  << std::endl;
    }


}


void test_richardson_extrapolation() {

    const double h_0 = 0.001;

    const double x_0 = 0.1;

    // f
    std::function<double(double)> f = [](double x) {
        return -1.0 / x;
    };

    // f'
    std::function<double(double)> f_ = [](double x) {
        return 1.0 / (x * x);
    };


    std::function<double(double)> f_approx = [&](double h) {
        return cfd::math::bnas::dfdx(x_0, h, f);
    };


    std::cout << "f'(x_0)[analytic] = f'(" << x_0 << ") = " << f_(x_0) << std::endl;

    std::cout << "~f'(x_0)[numeric] = ~f'(" << x_0 << ") = " << f_approx(h_0) << std::endl;

    std::cout << "~f'(x_0)[numeric,richardson] = ~f'(x_0)(" << x_0 << ") = "
              << cfd::math::bnas::richardson_extrapolation(h_0, f_approx, 10)
              << std::endl;
}

int main(int argc, const char *argv[]) {
    test_richardson_extrapolation();
    test_integrate_point();
    return 0;
}



//    u->forEach([&](auto &u_xy, auto &v_xy, int x, int y) {
//
//        auto x_n_plus_1_x = fluid::util::math::integration::rk3_x(f_x_q, x, y, h);
//        auto x_n_plus_1_y = fluid::util::math::integration::rk3_y(f_y_q, x, y, h);
//
//        double alpha = u_xy - x_n_plus_1_x;
//        double beta = v_xy - x_n_plus_1_y;
//
//        u_xy = fluid::util::math::interpolation::general::bi_apply_4th_order(-1, 0, 1, 2, x, y, f_x_q, kernel,
//                                                                                         alpha,
//                                                                                         beta,
//                                                                                         kernel_radius);
//
//        v_xy = fluid::util::math::interpolation::general::bi_apply_4th_order(-1, 0, 1, 2, x, y, f_y_q, kernel,
//                                                                                         alpha,
//                                                                                         beta,
//                                                                                         kernel_radius);
//    });
//
//    u->forEach([&](auto &u_xy, auto &v_xy, int x, int y) {
//
//        D(std::cout << "q(x,y) = q(" << x << "," << y << ") = " <<
//                    "(" << q->operator()(x, y).x() << "," << q->operator()(x, y).y() << ")" << std::endl;)
//
//        Vector2d x_n = Vector2d(x, y);
//        Vector2d q_x_n_plus_1;
//        Vector2d x_n_plus_1;
//
//#if FEATURE_ADVECT_WITH_RICHARDSON_EXTRAPOLATION == 0
//        x_n_plus_1.x() = fluid::util::math::integration::rk3_x(f_x_u, x, y, h);
//        x_n_plus_1.y() = fluid::util::math::integration::rk3_y(f_y_u, x, y, h);
//#else
//        // Richardson extrapolation (RE) yields optimized h which is much more precises. Doing it all the time might be costly.
//        // TODO Think of a process, running RE once to find stationary h
//        std::function<Vector2d(double)> integr_h = [&](double h) {
//            return fluid::util::math::integration::rk4(q, x_n, h);
//        };
//        x_n_plus_1 = cfd::math::bnas::richardson_extrapolation<Vector2d>(h, integr_h, 5);
//#endif
//
//        // round to 2 digits to cut off error
//#ifdef CUT_OFF_PRECISION
//        x_n_plus_1.x() = std::round(x_n_plus_1.x() * 1000.0) / 1000.0;
//        x_n_plus_1.y() = std::round(x_n_plus_1.y() * 1000.0) / 1000.0;
//#endif
//
//        // TODO Implement kernel radius specific scaling
//        // TODO Make vectorized / directional interpolation working!
//
//        double alpha = x_n.x() - x_n_plus_1.x();
//        double beta = x_n.y() - x_n_plus_1.y();
//
//        //q_x_n_plus_1.x() = fluid::util::math::interpolation::bridson(f_x, x, y, alpha, 0.0);
//        //q_x_n_plus_1.y() = fluid::util::math::interpolation::bridson(f_y, x, y, 0.0, beta);
//
//        //q_x_n_plus_1.x() = fluid::util::math::interpolation::bridson(f_x, x, y, alpha, beta);
//        //q_x_n_plus_1.y() = fluid::util::math::interpolation::bridson(f_y, x, y, alpha, beta);
//
////        q_x_n_plus_1.x() = fluid::util::math::interpolation::bilerp(x,  y,x + 1, y + 1, alpha, 0.0);
////        q_x_n_plus_1.y() = fluid::util::math::interpolation::bilerp(x,  y,x + 1, y + 1, 0.0, beta);
//
////        q_x_n_plus_1.x() = f_x(x, y);
////        q_x_n_plus_1.y() = f_y(x, y);
//
////        q_x_n_plus_1.x() = fluid::util::math::interpolation::general::bi_apply_4th_order(-1, 0, 1, 2, x, y, f_x_q, kernel,
////                                                                                         alpha,
////                                                                                         beta,
////                                                                                         kernel_radius);
////
////
////        q_x_n_plus_1.y() = fluid::util::math::interpolation::general::bi_apply_4th_order(-1, 0, 1, 2, x, y, f_y_q, kernel,
////                                                                                         alpha,
////                                                                                         beta,
////                                                                                         kernel_radius);
//         cell << q_x_n_plus_1;
//    });