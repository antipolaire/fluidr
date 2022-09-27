use vector2d::Vector2D;
//      const double DEFAULT_STEP_WIDTH = std::sqrt(std::numeric_limits<double>::epsilon());
//
//     /*
//     double richardson_extrapolation(double h, std::function<double(double)> N_1, int n = 2) {
//         using namespace std::placeholders;
//         std::function<double(std::function<double(double)>, double, int)>
//                 N_m = [](std::function<double(double)> N_m_minus_1, double h, int k) {
//             int pow_4_k = std::pow(4, k);
//             return (pow_4_k * N_m_minus_1(h) - N_m_minus_1(2.0 * h)) / (pow_4_k - 1);
//         };
//
//         std::function<double(double)> N_i = N_1;
//         for (int i = 2; i <= n; i++) {
//             N_i = std::bind(N_m, N_i, _1, i);
//         }
//         return N_i(h);
//     }
//     */

//
//     template<class T>
//     T richardson_extrapolation(double h, std::function<T(double)> N_1, int n = 2) {
//         using namespace std::placeholders;
//         std::function<T(std::function<T(double)>, double, int)>
//                 N_m = [](std::function<T(double)> N_m_minus_1, double h, int k) {
//             int pow_4_k = std::pow(4, k);
//             return (pow_4_k * N_m_minus_1(h) - N_m_minus_1(2.0 * h)) / (pow_4_k - 1);
//         };
//
//         std::function<T(double)> N_i = N_1;
//         for (int i = 2; i <= n; i++) {
//             N_i = std::bind(N_m, N_i, _1, i);
//         }
//         return N_i(h);
//     }

fn richardson_extrapolation<T, F>(h: f64, N_1: F, n: i32) -> T
where
    F: Fn(f64) -> T,
{
    let N_m = |N_m_minus_1: F, h: f64, k: i32| {
        let pow_4_k = 4_i32.pow(k as u32);
        (pow_4_k * N_m_minus_1(h) - N_m_minus_1(2.0 * h)) / (pow_4_k - 1)
    };

    let mut N_i = N_1;
    for i in 2..=n {
        N_i = |h| N_m(N_i, h, i);
    }
    N_i(h)
}

//
//     /**
//      * First order derivative of f at c with step width h
//      */
//     double dfdx(double x, double h, std::function<double(double)> f) {
//         return (f(x + h) - f(x - h)) / (2.0 * h);
//     }

fn dfdx(x: f64, h: f64, f: impl Fn(f64) -> f64) -> f64 {
    (f(x + h) - f(x - h)) / (2.0 * h)
}

//
//     /**
//      * Second order derivative of f at c with step width h
//      */
//     double d2fdx2(double x, double h, std::function<double(double)> f) {
//         return (-1.0 * f(x + h * 2.0) + 16.0 * f(x + h) - 30.0 * f(x) + 16.0 * f(x - h) - f(x - 2.0 * h)) /
//                (12 * h * h);
//         //return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
//     }

// Consider calling it with a step width of f64::EPSILON
fn d2fdx2(x: f64, h: f64, f: impl Fn(f64) -> f64) -> f64 {
    (-1.0 * f(x + h * 2.0) + 16.0 * f(x + h) - 30.0 * f(x) + 16.0 * f(x - h) - f(x - 2.0 * h)) / (12.0 * h * h)
}

//
//     double d2fdx2(double x, std::function<double(double)> f) {
//         return d2fdx2(x, std::sqrt(DEFAULT_STEP_WIDTH) * x, f);
//     }

fn d2fdx2_(x: f64, f: impl Fn(f64) -> f64) -> f64 {
    d2fdx2(x, f64::EPSILON * x, f)
}

//
//     double dfdx(double x, double y, double h, std::function<double(double, double)> f) {
//         return (f(x + h, y) - f(x - h, y)) / (2.0 * h);
//     }

fn dfdx_(x: f64, y: f64, h: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    (f(x + h, y) - f(x - h, y)) / (2.0 * h)
}

//
//     double dfdy(double x, double y, double h, std::function<double(double, double)> f) {
//         return (f(x, y + h) - f(x, y - h)) / (2.0 * h);
//     }

fn dfdy(x: f64, y: f64, h: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    (f(x, y + h) - f(x, y - h)) / (2.0 * h)
}
//
//     double dfdx(double x, double y, std::function<double(double, double)> f) {
//         return dfdx(x, y, DEFAULT_STEP_WIDTH * x, f);
//     }

fn dfdx__(x: f64, y: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    dfdx_(x, y, f64::EPSILON * x, f)
}

//
//     double dfdy(double x, double y, std::function<double(double, double)> f) {
//         return dfdy(x, y, DEFAULT_STEP_WIDTH * x, f);
//     }

fn dfdy_(x: f64, y: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    dfdy(x, y, f64::EPSILON * x, f)
}

//
//     double d2fdxx(double x, double y, double h, std::function<double(double, double)> f) {
//         return (f(x + h, y) - 2.0 * f(x, y) + f(x - h, y)) / (h * h);
//     }

fn d2fdxx(x: f64, y: f64, h: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    (f(x + h, y) - 2.0 * f(x, y) + f(x - h, y)) / (h * h)
}

//
//     double d2fdyy(double x, double y, double h, std::function<double(double, double)> f) {
//         return (f(x, y + h) - 2.0 * f(x, y) + f(x, y - h)) / (h * h);
//     }

fn d2fdyy(x: f64, y: f64, h: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    (f(x, y + h) - 2.0 * f(x, y) + f(x, y - h)) / (h * h)
}

//
//     double d2fdxy(double x, double y, double h, double k, std::function<double(double, double)> f) {
//         return (f(x + h, y + k) - f(x + h, y - k) - f(x - h, y + k) + f(x - h, y - k)) / (2.0 * h * k);
//     }

fn d2fdxy(x: f64, y: f64, h: f64, k: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    (f(x + h, y + k) - f(x + h, y - k) - f(x - h, y + k) + f(x - h, y - k)) / (2.0 * h * k)
}
//
//     double d2fdxy(double x, double y, std::function<double(double, double)> f) {
//         return d2fdxy(x, y, DEFAULT_STEP_WIDTH * x, DEFAULT_STEP_WIDTH * x, f);
//     }

fn d2fdxy_(x: f64, y: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    d2fdxy(x, y, f64::EPSILON * x, f64::EPSILON * x, f)
}

//
//     Eigen::Vector2d grad(double x, double y, std::function<double(double, double)> f) {
//         return Eigen::Vector2d(dfdx(x, y, f), dfdy(x, y, f));
//     }

fn grad(x: f64, y: f64, &f: impl Fn(f64, f64) -> f64) -> Vector2D<f64> {
    Vector2D::new(dfdx__(x, y, f), dfdy_(x, y, f))
}

//
//     double div(double x, double y, std::function<Eigen::Vector2d(double, double)> f) {
//         double h = DEFAULT_STEP_WIDTH * x;
//         return (f(x + h, y)(0) - f(x - h, y)(0)) / (2.0 * h) +
//                (f(x + h, y)(1) - f(x - h, y)(1)) / (2.0 * h);
//     }

fn div(x: f64, y: f64, f: impl Fn(f64, f64) -> Vector2D<f64>) -> f64 {
    let h = f64::EPSILON * x;
    (f(x + h, y)(0) - f(x - h, y)(0)) / (2.0 * h) + (f(x + h, y)(1) - f(x - h, y)(1)) / (2.0 * h)
}

//
//     /**
//      * Divergence of the gradient. A.k.a. ∆
//      */
//     double laplace(double x, double y, std::function<double(double, double)> f) {
//         using namespace std::placeholders;
//         return div(x, y, std::bind(grad, _1, _2, f));
//     }

fn laplace(x: f64, y: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    div(x, y, |x, y| grad(x, y, f))
}
