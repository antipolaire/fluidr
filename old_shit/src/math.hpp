#ifndef FLUID_SOLVER_RAW_MATH_HPP
#define FLUID_SOLVER_RAW_MATH_HPP

// #include "field2d.hpp"
#include "macros.hpp"
#include <Eigen/Core>

using namespace Eigen;
/**
 * Util functions collections some common math operations
 */
namespace fluid::util::math {
    /**
     * Articles:
     * http://www.ipol.im/pub/art/2011/g_lmii/
     */
    namespace interpolation {
        template<class TVal>
        TVal lerp(TVal x0, TVal x1, double alpha) {
            return (1.0 - alpha) * x0 + alpha * x1;
        }

        double cerp(double a, double b, double c, double d, double x) {
            double xsq = x*x;
            double xcu = xsq*x;

            double minV = std::min(a, std::min(b, std::min(c, d)));
            double maxV = std::max(a, std::max(b, std::max(c, d)));

            double t =
                    a*(0.0 - 0.5*x + 1.0*xsq - 0.5*xcu) +
                    b*(1.0 + 0.0*x - 2.5*xsq + 1.5*xcu) +
                    c*(0.0 + 0.5*x + 2.0*xsq - 1.5*xcu) +
                    d*(0.0 + 0.0*x - 0.5*xsq + 0.5*xcu);

            return std::min(std::max(t, minV), maxV);
        }

        template<class TVal>
        inline TVal bilerp(TVal v00, TVal v10, TVal v01, TVal v11, double s, double t) {
            const TVal a = lerp(v00, v01, t);
            const TVal b = lerp(v10, v11, t);
            return lerp(a, b, s);
        }

        constexpr double w_m1(double s) { return 0.5 * s * s - (1.0 / 6.0) * s * s * s - (1.0 / 3.0) * s; }

        constexpr double w_0(double s) { return 1.0 - s * s + 0.5 * s * s * s - 0.5 * s; }

        constexpr double w_1(double s) { return s + 0.5 * s * s - 0.5 * s * s * s; }

        constexpr double w_2(double s) { return (1.0 / 6.0) * s * s * s - (1.0 / 6.0) * s; }

        /**
         * estimate value 'q' at a fraction 's' between grid points 'x_{i}' and 'x_{i+1}
         */

        template<class TVal, class TOrd>
        TVal bridson(const std::function<TVal(TOrd, TOrd)> &f,
                     const TOrd x,
                     const TOrd y,
                     const double s_x = 0.5,
                     const double s_y = 0.5) {

            const auto p0_offset = -1;
            const auto p1_offset = 0;
            const auto p2_offset = 1;
            const auto p3_offset = 2;
            /* q_{j-1} */
            auto q_j_m1 = interpolation::w_m1(s_x) * f(x + p0_offset, y + p0_offset) +
                          interpolation::w_0(s_x) * f(x + p1_offset, y + p0_offset) +
                          interpolation::w_1(s_x) * f(x + p2_offset, y + p0_offset) +
                          interpolation::w_2(s_x) * f(x + p3_offset, y + p0_offset);
            /* q_{j} */
            auto q_j = interpolation::w_m1(s_x) * f(x + p0_offset, y + p1_offset) +
                       interpolation::w_0(s_x) * f(x + p1_offset, y + p1_offset) +
                       interpolation::w_1(s_x) * f(x + p2_offset, y + p1_offset) +
                       interpolation::w_2(s_x) * f(x + p3_offset, y + p1_offset);

            /* q_{j+1} */
            auto q_j_p1 = interpolation::w_m1(s_x) * f(x + p0_offset, y + p2_offset) +
                          interpolation::w_0(s_x) * f(x + p1_offset, y + p2_offset) +
                          interpolation::w_1(s_x) * f(x + p2_offset, y + p2_offset) +
                          interpolation::w_2(s_x) * f(x + p3_offset, y + p2_offset);
            /* q_{j+2} */
            auto q_j_p2 = interpolation::w_m1(s_x) * f(x + p0_offset, y + p3_offset) +
                          interpolation::w_0(s_x) * f(x + p1_offset, y + p3_offset) +
                          interpolation::w_1(s_x) * f(x + p2_offset, y + p3_offset) +
                          interpolation::w_2(s_x) * f(x + p3_offset, y + p3_offset);

            auto interp = interpolation::w_m1(s_y) * q_j_m1 +
                          interpolation::w_0(s_y) * q_j +
                          interpolation::w_1(s_y) * q_j_p1 +
                          interpolation::w_2(s_y) * q_j_p2;


            return interp;
            // return interp < 0.000000001 ? 0.0 : interp;
        }
/**
         * Kernel based implementation according to http://www.ipol.im/pub/art/2011/g_lmii/article.pdf
         */
        namespace general {


            template<class TVal>
            TVal scale_range(const TVal value, const TVal value_min, const TVal value_max, const TVal range_min,
                             const TVal range_max) {
                return (range_max - range_min) * ((value - value_min) / (value_max - value_min)) + range_min;
            }

            template<class TVal>
            std::function<TVal(TVal)>
            normalized_kernel(const std::function<TVal(TVal)> &kernel, const TVal kernel_radius,
                              const TVal value_min = -1.0, const TVal value_max = 1.0) {
                std::function<TVal(TVal)> wrapped = [=](TVal value) {
                    TVal value_norm = scale_range(value, value_min, value_max, -kernel_radius, kernel_radius);
                    TVal value_kernel = kernel(value_norm);
                    TVal value_denorm = scale_range(value_kernel, -kernel_radius, kernel_radius, value_min, value_max);

                    //std::cout << "\tDBG kernel: value orig         = " << value << std::endl;
                    //std::cout << "\tDBG kernel: value normalized   = " << value_norm << std::endl;
                    //std::cout << "\tDBG kernel: value kernel       = " << value_kernel << std::endl;
                    //std::cout << "\tDBG kernel: value denormalized = " << value_denorm << std::endl;

                    return value_denorm;
                };
                return wrapped;
            }


            /**
             * radius = 0.51f
             * @param x
             * @return
             */

            double kernel_nearest_neighbor(double x) {
                if (-0.5 <= x && x < 0.5)
                    return 1;
                else
                    return 0;
            }

            /**
             * radius = 1.0
             * @param x
             * @return
             */
            double kernel_bilinear(double x) {
                x = std::abs(x);

                if (x < 1.0)
                    return 1.0 - x;
                else
                    return 0.0;
            }

            /**
             * radius = 2.0
             * @param x
             * @return
             */
            double kernel_bicubic(double x) {
                // -0.5 yields third order accuracy (check: https://www.ipol.im/pub/art/2011/g_lmii/article.pdf)
                const double alpha = -0.5;
                x = std::abs(x);


                if (x < 2.0) {
                    if (x <= 1.0)
                        return ((alpha + 2.0) * x - (alpha + 3.0)) * x * x + 1.0;
                    else
                        return ((alpha * x - 5.0 * alpha) * x + 8.0 * alpha) * x - 4.0 * alpha;
                } else
                    return 0.0;
            }

            /**
             * radius = 2.0
             * @param x
             * @return
             */
            double kernel_lanczos2(double x) {
                if (-2 < x && x < 2) {
                    if (x != 0)
                        return sin(M_PI * x) * sin((M_PI / 2) * x) / ((M_PI * M_PI / 2) * x * x);
                    else
                        return 1;
                } else
                    return 0;
            }

            /**
             * radius = 3.0
             * @param x
             * @return
             */
            double kernel_lanczos3(double x) {
                if (-3 < x && x < 3) {
                    if (x != 0)
                        return sin(M_PI * x) * sin((M_PI / 3) * x) / ((M_PI * M_PI / 3) * x * x);
                    else
                        return 1;
                } else
                    return 0;
            }


            /**
             * radius = 4.0
             * @param x
             * @return
             */
            double kernel_lanczos4(double x) {
                if (-4 < x && x < 4) {
                    if (x != 0)
                        return sin(M_PI * x) * sin((M_PI / 4) * x) / ((M_PI * M_PI / 4) * x * x);
                    else
                        return 1;
                } else
                    return 0;
            }

            template<class TVal>
            TVal apply_2nd_order(const TVal fx1, const TVal fx2, const std::function<double(double)> &kernel,
                                 const double lambda,
                                 const double kernel_radius = 1.0) {
                auto _lambda = std::clamp(lambda, -1.0, 1.0);
                std::function<TVal(TVal)> _kernel = normalized_kernel(kernel, kernel_radius);
                return _kernel(-_lambda) * fx1 + _kernel(1 - _lambda) * fx2;
            }

            template<class TVal>
            TVal apply_4th_order(const TVal fx1, const TVal fx2, const TVal fx3, const TVal fx4,
                                 const std::function<TVal(TVal)> &kernel,
                                 const TVal lambda,
                                 const TVal kernel_radius = 2.0) {
                auto _lambda = std::clamp(lambda, -kernel_radius, kernel_radius);
                //auto _lambda = lambda;
                std::function<TVal(TVal)> _kernel = normalized_kernel(kernel, kernel_radius);
                return kernel(-1.0 - _lambda) * fx1 +
                       kernel(0.0 - _lambda) * fx2 +
                       kernel(1.0 - _lambda) * fx3 +
                       kernel(2.0 - _lambda) * fx4;
            }

            // TODO Put kernels along with according radius into a structure and pass that in as one
            // TODO There's a bug in here, find it!
            template<class TVal, class TOrd>
            TVal bi_apply_4th_order(const float offset1,
                                    const float offset2,
                                    const float offset3,
                                    const float offset4,
                                    const TOrd center_x,
                                    const TOrd center_y,
                                    const std::function<TVal(TOrd, TOrd)> &f,
                                    const std::function<TVal(TVal)> &kernel,
                                    const double _alpha = 0.5,
                                    const double _beta = 0.5,
                                    const double _kernel_radius = 2.0) {

                // double alpha = std::clamp(_alpha, -_kernel_radius, _kernel_radius);
                // double beta = std::clamp(_beta, -_kernel_radius, _kernel_radius);

                double alpha = _alpha;
                double beta = _beta;

                TVal fx1 = apply_4th_order(f(center_x + offset1, center_y + offset1),
                                           f(center_x + offset2, center_y + offset1),
                                           f(center_x + offset3, center_y + offset1),
                                           f(center_x + offset4, center_y + offset1), kernel, beta, _kernel_radius);

                TVal fx2 = apply_4th_order(f(center_x + offset1, center_y + offset2),
                                           f(center_x + offset2, center_y + offset2),
                                           f(center_x + offset3, center_y + offset2),
                                           f(center_x + offset4, center_y + offset2), kernel, beta, _kernel_radius);

                TVal fx3 = apply_4th_order(f(center_x + offset1, center_y + offset3),
                                           f(center_x + offset2, center_y + offset3),
                                           f(center_x + offset3, center_y + offset3),
                                           f(center_x + offset4, center_y + offset3), kernel, beta, _kernel_radius);

                TVal fx4 = apply_4th_order(f(center_x + offset1, center_y + offset4),
                                           f(center_x + offset2, center_y + offset4),
                                           f(center_x + offset3, center_y + offset4),
                                           f(center_x + offset4, center_y + offset4), kernel, beta, _kernel_radius);

                return apply_4th_order(fx1,
                                       fx2,
                                       fx3,
                                       fx4, kernel, alpha, _kernel_radius);
            }

        }


        double cubicInterpolate(double p[4], double x) {
            return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
                                                        x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
        }

        double bicubicInterpolate(double p[4][4], double x, double y) {
            double arr[4];
            arr[0] = cubicInterpolate(p[0], y);
            arr[1] = cubicInterpolate(p[1], y);
            arr[2] = cubicInterpolate(p[2], y);
            arr[3] = cubicInterpolate(p[3], y);
            return cubicInterpolate(arr, x);
        }
    }

    namespace integration {

        template<class TVal, class TOrd>
        TVal rk2(const std::function<TVal(TOrd, TOrd)> &f, const TOrd x, const TOrd y, double dt) {
            return x + dt * (x + 0.5 * dt * f(x, y));
        }

        /**
         * See "Runge-Kutta Methods with Minimum Error Bounds By Anthony Ralston" : https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
         */
//        Vector2d rk3(Field2D<Eigen::Vector2d> *&u, const Vector2d &x, double h) {
//            Vector2d k1 = h * u->operator()(x);
//            Vector2d k2 = h * u->operator()(static_cast<Vector2d>(x + 0.5 * h * k1));
//            Vector2d k3 = h * u->operator()(static_cast<Vector2d>(x + 0.75 * h * k2));
//            return x + h * ((2.0 / 9.0) * k1 + (1.0 / 3.0) * k2 + (4.0 / 9.0) * k3);
//        }
//
//        Vector2d rk3(Field2D<Eigen::Vector2d> *&u, int x, int y, double dt) {
//            return rk3(u, Vector2d(x, y), dt);
//        }


        template<class TVal, class TOrd>
        TVal rk3_x(const std::function<TVal(TOrd, TOrd)> &f, const TOrd x, TOrd y, TVal dt) {
            TVal k1 = f((x), (y));
            TVal k2 = f(std::nearbyint(x + 0.5 * dt * k1), (y));
            TVal k3 = f(std::nearbyint(x + (3.0 / 4.0) * dt * k2), (y));
            return x + dt * ((2.0 / 9.0) * k1 + (1.0 / 3.0) * k2 + (4.0 / 9.0) * k3);
        }

        template<class TVal, class TOrd>
        TVal rk3_y(const std::function<TVal(TOrd, TOrd)> &f, const TOrd x, TOrd y, TVal dt) {
            TVal k1 = f((x), (y));
            TVal k2 = f((x), std::nearbyint(y + 0.5 * dt * k1));
            TVal k3 = f((x), std::nearbyint(y + (3.0 / 4.0) * dt * k2));
            return y + dt * ((2.0 / 9.0) * k1 + (1.0 / 3.0) * k2 + (4.0 / 9.0) * k3);
        }

        Vector2d rk4(const std::function<double(int, int)> &f_u,
                     const std::function<double(int, int)> &f_v,
                     const int i, const int j, double h) {
            auto k1x = f_u(i, j);
            auto k1y = f_v(i, j);
            auto k2x = f_u(i + 0.5 * h * k1x, j + 0.5 * h * k1y);
            auto k2y = f_v(i + 0.5 * h * k1x, j + 0.5 * h * k1y);
            auto k3x = f_u(i + 0.5 * h * k2x, j + 0.5 * h * k2y);
            auto k3y = f_v(i + 0.5 * h * k2x, j + 0.5 * h * k2y);
            auto k4x = f_u(i + h * k3x, j + h * k3y);
            auto k4y = f_v(i + h * k3x, j + h * k3y);
            return Vector2d(i + (1.0 / 6.0) * h * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
                            j + (1.0 / 6.0) * h * (k1y + 2.0 * k2y + 2.0 * k3y + k4y));
        }

        Vector2d rk4(const std::function<Vector2d(int, int)> &f, const int i, int j, double h) {
            const std::function<double(int, int)> f_u = [&f](int x, int y) {
                return f(x, y).x();
            };
            const std::function<double(int, int)> f_v = [&f](int x, int y) {
                return f(x, y).y();
            };
            return rk4(f_u, f_v, i, j, h);
        }

// Compare http://www15.ovgu.de/ifme/l-numerik/mnmm-2-kapitel%206.pdf
        template<class TVal, class TOrd>
        TVal rk_iter(const std::function<TVal(TOrd, TOrd)> &f, TOrd x0, TOrd y0, TOrd x, TVal h) {
            int i, n = (int) ((x - x0) / h);
            auto y = y0;
            TVal k1, k2, k3, k4;
            for (i = 0; i <= n; i++) {
                k1 = h * f(x, y);
                k2 = h * f(x + 0.5 * h, y + 0.5 * k1);
                k3 = h * f(x + 0.5 * h, y + 0.5 * k2);
                k4 = h * f(x + h, y + k3);
                y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
                x0 = x0 + h;
            }
            return y;
        }

        template<class TVal, class TOrd>
        TVal rk4_x(const std::function<TVal(TOrd, TOrd)> &f, const TOrd x, TOrd y, double dt) {
            TVal k1 = f(x, y);
            TVal k2 = f(x + 0.5 * dt * k1, y);
            TVal k3 = f(x + 0.5 * dt * k2, y);
            TVal k4 = f(x + dt * k3, y);
            return x + (1.0 / 6.0) * dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }

        template<class TVal, class TOrd>
        TVal rk4_y(const std::function<TVal(TOrd, TOrd)> &f, const TOrd x, TOrd y, double dt) {
            TVal k1 = f(x, y);
            TVal k2 = f(x, y + 0.5 * dt * k1);
            TVal k3 = f(x, y + 0.5 * dt * k2);
            TVal k4 = f(x, y + dt * k3);
            return y + (1.0 / 6.0) * dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }

//        /**
//         * Method 5.15 of "Runge-Kutta Methods with Minimum Error Bounds By Anthony Ralston" : https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
//         */
//        Vector2d rk4_5_15(Field2D<Eigen::Vector2d> *&u, const Vector2d &x, double h) {
//            Vector2d k1 = u->operator()(x);
//            Vector2d k2 = u->operator()(static_cast<Vector2d>(x + Vector2d(0.4 * h, 0.4 * h * k1.y())));
//            Vector2d k3 = u->operator()(
//                    static_cast<Vector2d>(x + Vector2d(0.6 * h, -0.15 * h * k1.y() + 0.75 * h * k2.y())));
//            Vector2d k4 = u->operator()(static_cast<Vector2d>(x + Vector2d(h, (1.0 / 44.0) *
//                                                                              (19.0 * h * k1.y() -
//                                                                               15.0 * h * k2.y() +
//                                                                               40.0 * h * k3.y()))));
//            return x + h * (.17476028 * k1 - 0.55148066 * k2 + 1.20553560 * k3 + .17118478 * k4);
//        }
//
//        Vector2d rk4_5_15(Field2D<Eigen::Vector2d> *&u, int x, int y, double dt) {
//            return rk4_5_15(u, Vector2d(x, y), dt);
//        }
//
//        /**
//         * Method 5.12 of "Runge-Kutta Methods with Minimum Error Bounds By Anthony Ralston" : https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
//         */
//        Vector2d rk4_5_12(Field2D<Eigen::Vector2d> *&u, const Vector2d &x, double h) {
//            Vector2d k1 = u->operator()(x);
//            Vector2d k2 = u->operator()(static_cast<Vector2d>(x - Vector2d(0.4 * h, 0.4 * h * k1.y())));
//            Vector2d k3 = u->operator()(
//                    static_cast<Vector2d>(x -
//                                          Vector2d(0.45573725 * h,
//                                                   0.29697761 * h * k1.y() + 0.15875964 * h * k2.y())));
//            Vector2d k4 = u->operator()(static_cast<Vector2d>(x - Vector2d(h,
//                                                                           0.21810040 * h * k1.y() +
//                                                                           3.05096516 * h * k2.y() +
//                                                                           3.83286476 * h * k3.y())));
//            return x + h * (.17476028 * k1 - 0.55148066 * k2 + 1.20553560 * k3 + .17118478 * k4);
//        }
//
//        Vector2d rk4_5_12(Field2D<Eigen::Vector2d> *&u, int x, int y, double dt) {
//            return rk4_5_12(u, Vector2d(x, y), dt);
//        }
//
//        /**
//         * Runge-Kutta method to numerically integrate second order dlg's of the form f' = f(x,y).
//         * According to "E. Kamke, Differentialgleichungen - Lösungsmethoden und Lösungen, 1977, 142"
//         */
//        double rk4_kamke(const std::function<double(int, int)> &f, int x_i, int y_i, double h) {
//            auto k1 = h * f(x_i, y_i);
//            auto k2 = h * f(x_i + 0.5 * h, y_i + 0.5 * k1);
//            auto k3 = h * f(x_i + 0.5 * h, y_i + 0.5 * k2);
//            auto k4 = h * f(x_i + h, y_i + k3);
//            /* y_{i+1} - y_{i} = 1/6 * k_1 ...
//             * => y_{i+1} = y_{i} + 1/6 * k_1 ... */
//            auto y_diff_i_plus_1 = (1.0 / 6.0) * k1 + (1.0 / 3.0) * k2 + (1.0 / 3.0) * k3 + (1.0 / 6.0) * k4;
//            auto y_i_plus_1 = y_i + y_diff_i_plus_1;
//            return y_i_plus_1;
//        }
//
//        Vector2d rk4_kamke_2d(Field2D<Eigen::Vector2d> *&u, int x_i, int y_i, double h) {
//            Vector2d u_next;
//
//            u_next.x() = rk4_kamke([&](int x, int y) {
//                return u->operator()(x, y).x();
//            }, x_i, y_i, h);
//
//            u_next.y() = rk4_kamke([&](int x, int y) {
//                return u->operator()(x, y).y();
//            }, x_i, y_i, h);
//
//            return u_next;
//        }

    }

}
#endif //FLUID_SOLVER_RAW_MATH_HPP
