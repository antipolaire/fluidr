#ifndef FLUID_SOLVER_HPP
#define FLUID_SOLVER_HPP

#include <iostream>
#include <array>
#include <future>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <unistd.h>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
// #include <Eigen/CholmodSupport>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>

// TODO Either interpolation or RK4 is faulty or there is some weird numerical error. Maybe compare with code of some  other guys as https://github.com/AgentLee/FluidSim/

// Numeric optimization taken from https://d2f99xq7vri1nk.cloudfront.net/legacy_app_files/pdf/GDC03.pdf

// #define CUT_OFF_PRECISION

// #define DEBUG

#include "macros.hpp"
#include "bnas.h"

#include "field2d_v2.hpp"
#include "math.hpp"

using namespace Eigen;
using namespace std;

#define FEATURE_ADVECTION 1
// #define FEATURE_GRAVITY 1
#define FEATURE_PRESSURE 1
// #define FEATURE_VORTICITY 1
#define FEATURE_ADVECT_WITH_RICHARDSON_EXTRAPOLATION 0
// #define FEATURE_CENTER_OBSTACLE 1
// #define FEATURE_SOLID_WALLS 1

// #define INTEGRATION_RK4 1
// #define INTEGRATION_RK3 1
// #define INTEGRATION_STAM 1
#define INTEGRATION_EULER 1


#define IS_IN_RANGE(value, lower_limit, upper_limit) ((value)>(lower_limit) && (value)<(upper_limit))

class FluidSolver {
private:

    static const int CELL_TYPE_UNKNOWN = 0;
    static const int CELL_TYPE_FLUID = 1;
    static const int CELL_TYPE_SOLID = 2;
    static const int CELL_TYPE_EMPTY = 3;
    constexpr static const double CELL_SIZE = 0.5;
    constexpr static const double RK_H_DEFAULT = 2.0f;

public:
    double cos0 = 94.0f;
    double cos1 = 0.1f;
    double advectAmount = 1.0;
    double viscosity = 0.1;
    double diff = 0.0f;
    double dx = CELL_SIZE;
    double dt = 1.000;
    double rk_h = RK_H_DEFAULT;

    int width;
    int height;

    Field2DV2<double> *u;
    Field2DV2<double> *u_tmp;

    FluidSolver(int width, int height);

    virtual ~FluidSolver();

    void step();

    void fill(int style = 0);

    double getVelocityU(int x, int y) const;
    double getVelocityV(int x, int y) const;
    double getPressure(int x, int y) const;

    Field2DV2<double> *q;
private:

    const Vector2d SOLID_VELOCITY = Vector2d(0.0, 0.0);
    const Vector2d GRAVITY_CONSTANT = Vector2d(0.0, 0.91);

    int *cellType;

    std::function<void()> onRefreshData;

    using ivec2 = Eigen::Vector2i;
    using vec2 = Eigen::Vector2d;

    bool isCellType(int x, int y, int cellType);

    void advect(Field2DV2<double> *&u, double dt);

    void vorticity(Field2DV2<double> *&u, double dt, double vortCoeff);

    void update_pressure_gradient(double dt, double dx, double density, Field2DV2<double> *&u, VectorXd &p,
                                  int *&cellType);

    void init_negative_divergence(double dx, Field2DV2<double> *&u, int *&cellType,
                                  VectorXd &rhs);

    void init_pressure_equation_coefficients(double dt, double dx, double density, Field2DV2<double> *&u,
                                             int *&cellType, SparseMatrix<double> &A);

    void project(double dt, double dx, double viscosity, Field2DV2<double> *&u, int *&cellType);

    double curl(Field2DV2<double> *&u, int x, int y) const;

    void boundaries();

    double div(Field2DV2<double> *&u, double scale, int x, int y) const;

    void infuse(Field2DV2<double> *dst, Field2DV2<double> *src, double dt);

};

bool FluidSolver::isCellType(int x, int y, int type) {
    return this->cellType[XY2I(x,y,width)] == type;
}

/**
 * Idea:
 * Location in space of the grid point we're looking at is ̅x_G.
 * We want to find out the new value of q at that point, which we call q_G^{n+1}.
 * Go from ̅x_G backward in velocity ̅u_G with time-step ∆t to find the previous ̅x_P.
 * ̅x_P = ̅x_G − ∆t̅u_G.
 * q_G^{n+1} = interpolate (q^n, ̅x_G − ∆t̅u_G) = interpolate(q^n, ̅x_P)
 * @param u
 * @param dt
 * @param q
 */
void FluidSolver::advect(Field2DV2<double> *&u, double dt) {
    u->copyTo(u_tmp);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            if (!isCellType(i, j, CELL_TYPE_FLUID)) {
                continue;
            }

#ifdef INTEGRATION_RK4
            double k1x = this->u->u(i, j);
            double k1y = this->u->v(i, j);
            double _if = (double) i;
            double _jf = (double) j;
            double k2x = this->u->u(_if + 0.5f * rk_h * k1x, _jf + 0.5f * rk_h * k1y);
            double k2y = this->u->v(_if + 0.5f * rk_h * k1x, _jf + 0.5f * rk_h * k1y);
            double k3x = this->u->u(_if + 0.5f * rk_h * k2x, _jf + 0.5f * rk_h * k2y);
            double k3y = this->u->v(_if + 0.5f * rk_h * k2x, _jf + 0.5f * rk_h * k2y);
            double k4x = this->u->u(_if + rk_h * k3x, _jf + rk_h * k3y);
            double k4y = this->u->v(_if + rk_h * k3x, _jf + rk_h * k3y);
            double itrace = _if + (1.0f / 6.0f) * rk_h * (k1x + 2.0f * k2x + 2.0f * k3x + k4x);
            double jtrace = _jf + (1.0f / 6.0f) * rk_h * (k1y + 2.0f * k2y + 2.0f * k3y + k4y);

            // field is capable of interpolation
            this->u_new->u(i,j) = this->q->u(itrace,itrace);
            this->u_new->v(i,j) = this->q->v(itrace,jtrace);
#elif INTEGRATION_RK3
            double k1x = this->u->u(i, j);
            double k1y = this->u->v(i, j);
            double x = (double) i;
            double y = (double) j;
            double k2x = this->u->u(x + 0.5*dt*k1x, y + 0.5*dt*k1y);
            double k2y = this->u->v(x + 0.5*dt*k1x, y + 0.5*dt*k1y);
            double k3x = this->u->u(x + (3.0/4.0)*dt*k2x, y + (3.0/4.0)*dt*k2y);
            double k3y = this->u->v(x + (3.0/4.0)*dt*k2x, y + (3.0/4.0)*dt*k2y);
            u_tmp->u(i,j) = this->u->u(x + (2.0 / 9.0)*dt*k1x + (3.0 / 9.0)*dt*k2x + (4.0 / 9.0)*dt*k3x,
                                         y + (2.0 / 9.0)*dt*k1y + (3.0 / 9.0)*dt*k2y + (4.0 / 9.0)*dt*k3y);
            u_tmp->v(i,j) = this->u->v(x + (2.0 / 9.0)*dt*k1x + (3.0 / 9.0)*dt*k2x + (4.0 / 9.0)*dt*k3x,
                                         y + (2.0 / 9.0)*dt*k1y + (3.0 / 9.0)*dt*k2y + (4.0 / 9.0)*dt*k3y);
            // field is capable of interpolation
#elif INTEGRATION_EULER

            // integrate
            double x = i - dt * u->u(i, j);
            double y = j - dt * u->v(i, j);

            // interpolate
            if ( x < 0.5f ) {
                x = 0.5f;
            }
            if ( x > width + 0.5f ) {
                x = width + 0.5f;
            }
            int i0 = (int) x;
            int i1 = i0 + 1;
            if ( y < 0.5f ) {
                y = 0.5f;
            }
            if ( y > height + 0.5f ) {
                y = height + 0.5f;
            }
            int j0 = (int) y;
            int j1 = j0 + 1;
            double s1 = x - i0;
            double s0 = 1 - s1;
            double t1 = y - j0;
            double t0 = 1 - t1;

            u_tmp->u(i,j) = s0 * (t0 * q->u(i0, j0) + t1 * q->u(i0, j1)) +
                   s1 * (t0 * q->u(i1, j0) + t1 * q->u(i1, j1));

            u_tmp->v(i,j) = s0 * (t0 * q->v(i0, j0) + t1 * q->v(i0, j1)) +
                   s1 * (t0 * q->v(i1, j0) + t1 * q->v(i1, j1));
#endif
        }
    }
    u_tmp->copyTo(u);
}

void FluidSolver::vorticity(Field2DV2<double> *&u, double dt, double vortCoeff){
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            if (!isCellType(x, y, CELL_TYPE_FLUID)) {
                continue;
            }

            double dx = abs(curl(u, x, y - 1)) - abs(curl(u, x, y + 1));
            double dy = abs(curl(u, x + 1, y)) - abs(curl(u, x - 1, y));
            double len =sqrt(sqrt(dx)+sqrt(dy));
            dx = vortCoeff / len * dx;
            dy = vortCoeff / len * dy;
            u->u(x,y) = u->u(x,y) + dt * curl(u, x,y) * dx;
            u->v(x,y) = u->v(x,y) + dt * curl(u, x,y) * dy;
        }
    }
}

double FluidSolver::curl(Field2DV2<double> *&u, int x, int y) const {
    return u->u(x, y + 1) - u->u(x, y - 1) + u->v(x - 1, y) - u->v(x + 1, y);
}

void FluidSolver::project(double dt, double dx, double viscosity, Field2DV2<double> *&u, int *&cellType) {

    int long square_size = u->getSize();

    // rhs is the right hand side "b" in "Ap = b"
    VectorXd rhs(square_size);
    VectorXd p(square_size);

    // 1. Calculate negative divergence b (the right hand side) with modifications at solid wall boundaries (Figure 5.3)
    init_negative_divergence(dx, u, cellType, rhs);

    SparseMatrix<double> A(square_size, square_size);

    A.reserve(VectorXd::Constant( square_size,6));

    init_pressure_equation_coefficients(dt, dx, viscosity, u, cellType, A);

    D(std::cout << "A = " << std::endl << A << std::endl;)

    Eigen::ConjugateGradient<SparseMatrix<double>> solver(A);

    p = solver.solve(rhs);

    // 5. Compute new velocities ̅u^{n+1} according to pressure-gradient update to ̅u
    update_pressure_gradient(dt, dx, viscosity, u, p, cellType);

}

void FluidSolver::init_pressure_equation_coefficients(double dt, double dx, double density, Field2DV2<double> *&u,
                                                      int *&cellType, SparseMatrix<double> &A) {
    double scale = dt / (density * dx * dx);

    int i = 0,
     diag_inc = 0,
     right_inc = 0,
     bottom_inc = 0,
     top_inc = 0,
     left_inc = 0;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (!isCellType(x, y, CELL_TYPE_FLUID)) {
                continue;
            }

            i = XY2I(x,y,width);

            diag_inc = 1;
            right_inc = 0;
            bottom_inc = -1;
            top_inc = 0;
            left_inc = -1;
            if (x < width - 1) {
                diag_inc++;
                right_inc++;
                if(left_inc == -1){
                    left_inc = 0;
                }
                left_inc++;
            } else{
                left_inc = 0.0;
            }

            if (y < height - 1) {
                diag_inc++;
                top_inc++;
                if(bottom_inc == -1){
                    bottom_inc = 0;
                }
                bottom_inc++;
            } else{
                bottom_inc = 0;
            }

            if (diag_inc > 0) {
                A.coeffRef(i, i) = scale * diag_inc;
            }
//            if (right_inc > 0) {
//                A.coeffRef(i, i+1) = -scale * right_inc;
//            }
//            if (top_inc > 0) {
//                A.coeffRef(i, XY2I(x,y-1,width)) = scale * top_inc;
//            }
//            if (left_inc > -1) {
//                A.coeffRef(i, i-1) = -scale * left_inc;
//            }
//            if (bottom_inc > -1) {
//                A.coeffRef(i, XY2I(x,y+1,width)) = scale * bottom_inc;
//            }

        }
    }
    // A.makeCompressed(); // save space but slightly slower!
    return;
}

/**
 * Same as negative_divergence but with modification to respect solid walls
 * @param dt
 * @param dx
 * @param u
 * @param cellType
 * @param rhs
 */
void FluidSolver::init_negative_divergence(double dx, Field2DV2<double> *&u, int *&cellType,
                                           VectorXd &rhs) {

    double scale = 1.0 / dx;
    double neg_divergence_value;

    D(std::cout << "init negative divergence -∇̅̅u(length=" << rhs.size() << ")" << std::endl;)
    D(std::cout << std::fixed);
    D(std::cout << std::setprecision(3));
    // #pragma omp parallel for
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            neg_divergence_value = 0.0;
            if (isCellType( x, y, CELL_TYPE_FLUID)) {

                neg_divergence_value = div(u, scale, x, y) * -1.0;
                // Following is required for correct boundary conditions
                if (isCellType( x - 1, y, CELL_TYPE_SOLID)) {
                    neg_divergence_value = div(u, scale, x-1, y) * -1.0;
                }
                if (isCellType( x + 1, y, CELL_TYPE_SOLID)) {
                    neg_divergence_value = div(u, scale, x+1, y) * -1.0;
                }

                if (isCellType( x, y - 1, CELL_TYPE_SOLID)) {
                    neg_divergence_value = div(u, scale, x, y-1) * -1.0;
                }
                if (isCellType( x, y + 1, CELL_TYPE_SOLID)) {
                    neg_divergence_value = div(u, scale, x, y+1) * -1.0;
                }
            }

            D(std::cout << "rhs_xy("<< x << ","<< y << ") = rhs_index(" << XY2I(x,y,width) << ") = "<< neg_divergence_value << std::endl;)
            rhs(XY2I(x,y,width)) = neg_divergence_value;
        }
    }
}

double FluidSolver::div(Field2DV2<double> *&u, double scale, int x, int y) const {
    return ((u->u(x + 1, y) - u->u(x, y)) * scale + (u->v(x, y + 1) - u->v(x, y)) * scale);
}

void FluidSolver::update_pressure_gradient(double dt, double dx, double density, Field2DV2<double> *&u, VectorXd &p,
                                           int *&cellType) {
    const double scale = dt / (density * dx);
    D(std::cout << std::fixed);
    D(std::cout << std::setprecision(3));

#ifdef DEBUG
    int i = 0;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            i = XY2I(y, x, width);
            D(std::cout << "p(" << x << "," << y << ") = p("<< i <<") = " << p[i] << " | ");
        }
        D(std::cout << std::endl);
    }
#endif

    int idx=0;
    for ( int x = 0; x < width; x++ ) {
        for ( int y = 0; y < height; y++,idx++) {
            if (!isCellType(x, y, CELL_TYPE_FLUID)) {
                continue;
            }
            D(std::cout << "∇p(" << x << "," << y << ") : " << "(" << u->u(x,y) << "," << u->v(x,y) << ")" << std::endl);
            // The following approach is slower but better for complex boundary conditions (e.g. obstacles):
            if ((isCellType( x, y, CELL_TYPE_FLUID) || isCellType( x - 1, y, CELL_TYPE_FLUID))) {
                if (isCellType( x, y, CELL_TYPE_SOLID) || isCellType( x - 1, y, CELL_TYPE_SOLID)) {
                    // SOLID_VELOCITY
                    u->u(x,y) = 0.0;
                    u->v(x,y) = 0.0;
                } else {
                    u->u(x,y) -= scale * (p[XY2I(x, y,width)] - p[XY2I(x - 1, y, width)]);
                }
            } else {
                cellType[XY2I(x,y,width)] = CELL_TYPE_UNKNOWN;
            }

            if ((isCellType( x, y, CELL_TYPE_FLUID) || isCellType( x, y - 1, CELL_TYPE_FLUID))) {
                if (isCellType( x, y, CELL_TYPE_SOLID) || isCellType( x, y - 1, CELL_TYPE_SOLID)) {
                    // SOLID_VELOCITY
                   u->u(x,y) = 0.0;
                   u->v(x,y) = 0.0;
                } else {
                    u->v(x,y) -= scale * (p[XY2I(x, y,width)] - p[XY2I(x , y - 1, width)]);
                }
            } else {
                cellType[XY2I(x,y,width)] = CELL_TYPE_UNKNOWN;
            }

        }
    }

}

FluidSolver::FluidSolver(int width, int height) {
    this->width = width;
    this->height = height;
    this->dx = 1.0 / std::min(width, width);
    this->rk_h = this->dt;
    this->u = new Field2DV2<double>(width, height);
    this->u_tmp = new Field2DV2<double>(width, height);
    this->q = new Field2DV2<double>(width, height);

    this->cellType = new int[width * height];

    for ( int x = 0; x < width; x++ ) {
        for ( int y = 0; y < height; y++ ) {
            if (x == 0 || y == 0 || x == width - 1 || y == height - 1) {
#ifdef FEATURE_SOLID_WALLS
                cellType[XY2I(x,y,width)] = CELL_TYPE_SOLID;
#endif
            } else {
                cellType[XY2I(x,y,width)] = CELL_TYPE_FLUID;
            }
        }
    }
}

FluidSolver::~FluidSolver() {
    delete u;
    delete q;
    delete cellType;
}


void FluidSolver::fill(int style){
    int center_y = height / 2;
    switch (style) {
        case 0:
            for(int i=0;i<width*height;i++){
                q->u(i) = q->v(i) = 0.0;
            }
            for ( int y = 10; y < height-10; y++) {
                if(std::cos(((double)y)*cos0)*M_PI>cos1){
                    q->u(3,y) = advectAmount;
                }
            }
            break;
        case 1:
            for(int i=0;i<width*height;i++){
                q->u(i) = q->v(i) = 0.0;
            }

            q->u(3, center_y - 1) = advectAmount;
            q->u(3, center_y) = advectAmount;
            q->u(3, center_y + 1) = advectAmount;
            break;
        case 2:
            for(int i=0;i<width*height;i++){
                q->u(i) = q->v(i) = 0.0;
            }
            /* put some obstacle */
            for (int wall_idx = -30; wall_idx < 30; wall_idx++) {
                cellType[XY2I((width / 2 - 1) + 1, (center_y - 1) + wall_idx, width)] = CELL_TYPE_SOLID;
            }
            break;
        case 3:
            for(int i=0;i<width*height;i++){
                q->u(i) = q->v(i) = 0.0;
            }
            q->u((width / 2)+4, center_y) = advectAmount;
            q->u((width / 2)-4, center_y) = -advectAmount;

            q->v((width / 2), center_y + 4) = -advectAmount;
            q->v((width / 2), center_y - 4) = advectAmount;
            break;
    }
    infuse(u, q, dt);
}

void FluidSolver::step() {

#ifdef FEATURE_ADVECTION

#ifdef DEBUG
    std::clock_t c_start = std::clock();
#endif
    advect(u, dt);
#endif
#ifdef FEATURE_VORTICITY
    vorticity(u,dt,100.0);
#endif
#ifdef FEATURE_GRAVITY
    u->forEach([&](auto &cell, int x, int y) {
        cell += dt * GRAVITY_CONSTANT;
    });
#endif

#ifdef FEATURE_PRESSURE
    project(dt, dx, viscosity, u, cellType);
#endif
    boundaries();
#ifdef DEBUG
    std::clock_t c_end = std::clock();
    long time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "cpu time: " << time_elapsed_ms << " ms" << std::endl;
#endif
    // u->copyTo(q);
}

void FluidSolver::boundaries() {
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            if (!isCellType(x, y, CELL_TYPE_FLUID)) {
                this->u->u(x,y) = 0.0f;
                this->u->v(x,y) = 0.0f;
            }
        }
    };
}

double FluidSolver::getVelocityU(int x, int y) const {
    return u->u(x,y);
}

double FluidSolver::getVelocityV(int x, int y) const{
    return u->v(x,y);
}

double FluidSolver::getPressure(int x, int y) const {
    // TODO Hand out pressure later!
    return 0.0;
}

void FluidSolver::infuse(Field2DV2<double> *dst, Field2DV2<double> *src, double dt) {
    for(int i=0;i<width*height;i++){
        dst->u(i) += src->u(i);
        dst->v(i) += src->v(i);
    }
}


#endif //FLUID_SOLVER_RAW_FIELD2D_HPP
