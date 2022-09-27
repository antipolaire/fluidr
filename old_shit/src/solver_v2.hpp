#include <cfloat>
#include "../src/solver.hpp"

using namespace std;

#ifndef SOLVERH
/**
 * Fluid simulation / solver based on the jos stam's book "The Art of Fluid Animation". The unmodified core
 * implementation of this solver is intended to be used as reference to compare results of the robert bridson based
 * solver.
 */
class FluidSolverJosStam {

public:

    FluidSolverJosStam(int width, int height);

    virtual ~FluidSolverJosStam();

    void printVelocity();

    float getVelocityU(int x, int y);
    float getVelocityV(int x, int y);
    float getPressure(int x, int y);

    void step();

    void fill(int style = 0);

    float cos0 = 100.0f;
    float cos1 = 0.45f;
    float advectAmount = 10.0;
    float viscosity = 0.00001;
    float density = 0.0f;
    float dt = 0.05;

private:

    int width = 0;
    int height = 0;

    float * u;
    float * v;
    float * u0;
    float * v0;
    float * dens;
    float * dens_prev;

    void add_source(int N, float *&x, float *&s, float dt);

    void set_bnd(int N, int b, float *&x);

    void lin_solve(int N, int b, float *&x, float *&x0, float a, float c);

    void dens_step(int N, float *&x, float *&x0, float *&u, float *&v, float diff, float dt);

    void vel_step(int N, float *&u, float *&v, float *&u0, float *&v0, float viscosity, float dt);

    void project(int N, float *&u, float *&v, float *&p, float *&div);

    void advect(int N, int b, float *&d, float *&d0, float *&u, float *&v, float dt);

    void diffuse(int N, int b, float *&x, float *&x0, float viscosity, float dt);

};


#define SOLVERH

#define IX(i, j) ((i)+(N+1)*(j))
#define SWAP(x0, x) {float * tmp=x0;x0=x;x=tmp;}

void FluidSolverJosStam::add_source(int N, float *&x, float *&s, float dt) {
    int i, size = (N + 2) * (N + 2);
    for ( i = 0; i < size; i++ ) {
        x[i] += dt * s[i];
    }
}
void FluidSolverJosStam::set_bnd(int N, int b, float *&x) {
    int i;
    for ( i = 1; i <= N; i++ ) {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void FluidSolverJosStam::lin_solve(int N, int b, float *&x, float *&x0, float a, float c) {
    int i, j, k;

    for ( k = 0; k < 2; k++ ) {
        for ( i = 1; i <= N; i++ ) {
            for ( j = 1; j <= N; j++ ) {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
            }
        }
        set_bnd(N, b, x);
    }
}

void FluidSolverJosStam::diffuse(int N, int b, float *&x, float *&x0, float viscosity, float dt) {
    float a = dt * viscosity * N * N;
    lin_solve(N, b, x, x0, a, 1 + 4 * a);
}

void FluidSolverJosStam::advect(int N, int b, float *&d, float *&d0, float *&u, float *&v, float dt) {
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * N;
    for ( i = 1; i <= N; i++ ) {
        for ( j = 1; j <= N; j++ ) {
            x = i - dt0 * u[IX(i, j)];
            y = j - dt0 * v[IX(i, j)];
            if ( x < 0.5f ) {
                x = 0.5f;
            }
            if ( x > N + 0.5f ) {
                x = N + 0.5f;
            }
            i0 = (int) x;
            i1 = i0 + 1;
            if ( y < 0.5f ) {
                y = 0.5f;
            }
            if ( y > N + 0.5f ) {
                y = N + 0.5f;
            }
            j0 = (int) y;
            j1 = j0 + 1;
            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;
            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                          s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(N, b, d);
}

void FluidSolverJosStam::project(int N, float *&u, float *&v, float *&p, float *&div) {
    int i, j;

    for ( i = 1; i <= N; i++ ) {
        for ( j = 1; j <= N; j++ ) {
            div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(N, 0, div);
    set_bnd(N, 0, p);

    lin_solve(N, 0, p, div, 1, 4);

    for ( i = 1; i <= N; i++ ) {
        for ( j = 1; j <= N; j++ ) {
            u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
            v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
        }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);
}

void FluidSolverJosStam::dens_step(int N, float *&x, float *&x0, float *&u, float *&v, float diff, float dt) {
    add_source(N, x, x0, dt);
    SWAP (x0, x);
    diffuse(N, 0, x, x0, diff, dt);
    SWAP (x0, x);
    advect(N, 0, x, x0, u, v, dt);
}

void FluidSolverJosStam::vel_step(int N, float *&u, float *&v, float *&u0, float *&v0, float viscosity, float dt) {
    add_source(N, u, u0, dt);
    add_source(N, v, v0, dt);
    SWAP (u0, u);
    diffuse(N, 1, u, u0, viscosity, dt);
    SWAP (v0, v);
    diffuse(N, 2, v, v0, viscosity, dt);
    project(N, u, v, u0, v0);
    SWAP (u0, u);
    SWAP (v0, v);
    advect(N, 1, u, u0, u0, v0, dt);
    advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}

void FluidSolverJosStam::printVelocity(){
    std::cout << "u:" << std::endl;
    for (int i = 1; i <= this->width; ++i) {
        for (int j = 1; j <= this->height; ++j) {
            std::cout << "(" << getVelocityU(i,j) << "," << getVelocityV(i,j) << ")";
        }
        std::cout << std::endl;
    }
}

void FluidSolverJosStam::fill(int direction){
    int N = width;
    for ( int y = 10; y < height-10; y++) {
        if(std::cos(((float)y)*cos0)>cos1){
            this->u0[IX(3, y)] = advectAmount;
        }
    }
}

void FluidSolverJosStam::step(){
    vel_step(width, u, v, u0, v0, viscosity, dt);
    dens_step (width, dens, dens_prev, u, v, density, dt );
}

float FluidSolverJosStam::getVelocityU(int x, int y){
    int N = width;
    return u[IX(x,y)];
}
float FluidSolverJosStam::getVelocityV(int x, int y){
    int N = width;
    return v[IX(x,y)];
}

float FluidSolverJosStam::getPressure(int x, int y){
    int N = width;
    return dens[IX(x,y)];
}

FluidSolverJosStam::FluidSolverJosStam(int width, int height) {
    this->width = width;
    this->height = height;
    int size = (width+2)*(height+2);
    u			= (float *) malloc ( size*sizeof(float) );
    v			= (float *) malloc ( size*sizeof(float) );
    u0		    = (float *) malloc ( size*sizeof(float) );
    v0		    = (float *) malloc ( size*sizeof(float) );
    dens		= (float *) malloc ( size*sizeof(float) );
    dens_prev	= (float *) malloc ( size*sizeof(float) );
    int N = width;
}

FluidSolverJosStam::~FluidSolverJosStam() {
    if ( u ) free ( u );
    if ( v ) free ( v );
    if ( u0 ) free ( u0 );
    if ( v0 ) free ( v0 );
    if ( dens ) free ( dens );
    if ( dens_prev ) free ( dens_prev );
}

#endif