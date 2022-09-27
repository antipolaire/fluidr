//
// Created by Mona Lisa Overdrive
//

#ifndef CFD_CMAES_GL_BNASUTIL_H
#define CFD_CMAES_GL_BNASUTIL_H

#include <iostream>
#include <array>
#include <future>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <unistd.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>

/**
 * Basic Numerical Analysis Subprograms
 */
namespace cfd::math::bnas {

    const double DEFAULT_STEP_WIDTH = std::sqrt(std::numeric_limits<double>::epsilon());

    /*
    double richardson_extrapolation(double h, std::function<double(double)> N_1, int n = 2) {
        using namespace std::placeholders;
        std::function<double(std::function<double(double)>, double, int)>
                N_m = [](std::function<double(double)> N_m_minus_1, double h, int k) {
            int pow_4_k = std::pow(4, k);
            return (pow_4_k * N_m_minus_1(h) - N_m_minus_1(2.0 * h)) / (pow_4_k - 1);
        };

        std::function<double(double)> N_i = N_1;
        for (int i = 2; i <= n; i++) {
            N_i = std::bind(N_m, N_i, _1, i);
        }
        return N_i(h);
    }
    */

    template<class T>
    T richardson_extrapolation(double h, std::function<T(double)> N_1, int n = 2) {
        using namespace std::placeholders;
        std::function<T(std::function<T(double)>, double, int)>
                N_m = [](std::function<T(double)> N_m_minus_1, double h, int k) {
            int pow_4_k = std::pow(4, k);
            return (pow_4_k * N_m_minus_1(h) - N_m_minus_1(2.0 * h)) / (pow_4_k - 1);
        };

        std::function<T(double)> N_i = N_1;
        for (int i = 2; i <= n; i++) {
            N_i = std::bind(N_m, N_i, _1, i);
        }
        return N_i(h);
    }

    /**
     * First order derivative of f at c with step width h
     */
    double dfdx(double x, double h, std::function<double(double)> f) {
        return (f(x + h) - f(x - h)) / (2.0 * h);
    }

    /**
     * Second order derivative of f at c with step width h
     */
    double d2fdx2(double x, double h, std::function<double(double)> f) {
        return (-1.0 * f(x + h * 2.0) + 16.0 * f(x + h) - 30.0 * f(x) + 16.0 * f(x - h) - f(x - 2.0 * h)) /
               (12 * h * h);
        //return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
    }

    double d2fdx2(double x, std::function<double(double)> f) {
        return d2fdx2(x, std::sqrt(DEFAULT_STEP_WIDTH) * x, f);
    }

    double dfdx(double x, double y, double h, std::function<double(double, double)> f) {
        return (f(x + h, y) - f(x - h, y)) / (2.0 * h);
    }

    double dfdy(double x, double y, double h, std::function<double(double, double)> f) {
        return (f(x, y + h) - f(x, y - h)) / (2.0 * h);
    }

    double dfdx(double x, double y, std::function<double(double, double)> f) {
        return dfdx(x, y, DEFAULT_STEP_WIDTH * x, f);
    }

    double dfdy(double x, double y, std::function<double(double, double)> f) {
        return dfdy(x, y, DEFAULT_STEP_WIDTH * x, f);
    }

    double d2fdxx(double x, double y, double h, std::function<double(double, double)> f) {
        return (f(x + h, y) - 2.0 * f(x, y) + f(x - h, y)) / (h * h);
    }

    double d2fdyy(double x, double y, double h, std::function<double(double, double)> f) {
        return (f(x, y + h) - 2.0 * f(x, y) + f(x, y - h)) / (h * h);
    }

    double d2fdxy(double x, double y, double h, double k, std::function<double(double, double)> f) {
        return (f(x + h, y + k) - f(x + h, y - k) - f(x - h, y + k) + f(x - h, y - k)) / (2.0 * h * k);
    }

    double d2fdxy(double x, double y, std::function<double(double, double)> f) {
        return d2fdxy(x, y, DEFAULT_STEP_WIDTH * x, DEFAULT_STEP_WIDTH * x, f);
    }

    Eigen::Vector2d grad(double x, double y, std::function<double(double, double)> f) {
        return Eigen::Vector2d(dfdx(x, y, f), dfdy(x, y, f));
    }

    double div(double x, double y, std::function<Eigen::Vector2d(double, double)> f) {
        double h = DEFAULT_STEP_WIDTH * x;
        return (f(x + h, y)(0) - f(x - h, y)(0)) / (2.0 * h) +
               (f(x + h, y)(1) - f(x - h, y)(1)) / (2.0 * h);
    }

    /**
     * Divergence of the gradient. A.k.a. ∆
     */
    double laplace(double x, double y, std::function<double(double, double)> f) {
        using namespace std::placeholders;
        return div(x, y, std::bind(grad, _1, _2, f));
    }

    /**
     * D̅x/Dt = ∂̅̅x/∂t + (̅u⋅∇)̅x.
     * =>
     * D̅x/Dt = ∂̅̅x/∂t + (̅u⋅∇)̅x = ∂̅̅x/∂t + u₁*∂̅̅x/∂x₁ + u₂*∂̅̅x/∂x₂ + u₃*∂̅̅x/∂x₃.
     *
     * (x,y,z) -> [x,y,z]
     *
     */


    void test() {
        std::function<Eigen::Vector2d(double, double, int)> u;
        double x1, x2;
        int t;
        double u1, u2;

        Eigen::Vector2d mat_deriv =
                (u(x1, x2, t + 1) - u(x1, x2, t - 1)) / 2.0 +
                u1 * ((u(x1 + 1, x2, t) - u(x1 - 1, x2, t)) / 2.0) +
                u2 * ((u(x1, x2 + 1, t) - u(x1, x2 - 1, t)) / 2.0);
    }

    /*
    Eigen::Vector2d
    material_derivative(double x, double y, Eigen::Vector2d u, std::function<double(double, double)> f) {
        return dfdx(x, y, DEFAULT_STEP_WIDTH, f) + (u * grad(x, y, f));
    }
    */

    void curl() {

    }
}


#endif //CFD_CMAES_GL_BNASUTIL_H
