// /**
//  * Fluid simulation / solver based on the jos stam's book "The Art of Fluid Animation". The unmodified core
//  * implementation of this solver is intended to be used as reference to compare results of the robert bridson based
//  * solver.
//  */
// class FluidSolverJosStam {

use vector2d::Vector2D;

// #define IX!(i, j) ((i)+(N+1)*(j))

#[macro_export]
macro_rules! IX {
    ($i:expr, $j:expr) => {
        ($i) + ($N + 1) * ($j)
    };
}

// #define SWAP(x0, x) {float * tmp=x0;x0=x;x=tmp;}

#[macro_export]
macro_rules! SWAP {
    ($x0:expr, $x:expr) => {
        let tmp = $x0;
        $x0 = $x;
        $x = tmp;
    };
}


trait Solver {
    fn new() -> Self;
    fn step(&mut self, dt: f32);
    fn add_density(&mut self, x: usize, y: usize, amount: f32);
    fn add_velocity(&mut self, x: usize, y: usize, amount_x: f32, amount_y: f32);
    fn get_density(&self, x: usize, y: usize) -> f32;
    fn get_velocity(&self, x: usize, y: usize) -> (f32, f32);

    fn print_velocity(&self);
    fn get_velocity_u(&self, x: usize, y: usize) -> f32;
    fn get_velocity_v(&self, x: usize, y: usize) -> f32;
    fn get_pressure(&self, x: usize, y: usize) -> f32;
    fn step_(&mut self);

    fn fill(&mut self, style: i32);

    //
//     float cos0 = 100.0f;
//     float cos1 = 0.45f;
//     float advectAmount = 10.0;
//     float viscosity = 0.00001;
//     float density = 0.0f;
//     float dt = 0.05;

    const cos0: f32 = 100.0;
    const cos1: f32 = 0.45;
    const advectAmount: f32 = 10.0;
    const viscosity: f32 = 0.00001;
    const density: f32 = 0.0;
    const dt: f32 = 0.05;
}

//
// public:
//
//     FluidSolverJosStam(int width, int height);

//
//     virtual ~FluidSolverJosStam();
//
//     void printVelocity();
//
//     float getVelocityU(int x, int y);
//     float getVelocityV(int x, int y);
//     float getPressure(int x, int y);
//
//     void step();
//
//     void fill(int style = 0);
//
//     float cos0 = 100.0f;
//     float cos1 = 0.45f;
//     float advectAmount = 10.0;
//     float viscosity = 0.00001;
//     float density = 0.0f;
//     float dt = 0.05;
//
// private:
//
//     int width = 0;
//     int height = 0;
//
//     float * u;
//     float * v;
//     float * u0;
//     float * v0;
//     float * dens;
//     float * dens_prev;
//
//     void add_source(int N, float *&x, float *&s, float dt);
//
//     void set_bnd(int N, int b, float *&x);
//
//     void lin_solve(int N, int b, float *&x, float *&x0, float a, float c);
//
//     void dens_step(int N, float *&x, float *&x0, float *&u, float *&v, float diff, float dt);
//
//     void vel_step(int N, float *&u, float *&v, float *&u0, float *&v0, float viscosity, float dt);
//
//     void project(int N, float *&u, float *&v, float *&p, float *&div);
//
//     void advect(int N, int b, float *&d, float *&d0, float *&u, float *&v, float dt);
//
//     void diffuse(int N, int b, float *&x, float *&x0, float viscosity, float dt);
//
// };

impl dyn Solver {
    fn add_density(&mut self, x: usize, y: usize, amount: f32) {
        self.dens[x + y * self.width] += amount;
    }

    fn add_velocity(&mut self, x: usize, y: usize, amount_x: f32, amount_y: f32) {
        let index = x + y * self.width;
        self.u[index] += amount_x;
        self.v[index] += amount_y;
    }

    fn get_density(&self, x: usize, y: usize) -> f32 {
        self.dens[x + y * self.width]
    }

    fn get_velocity(&self, x: usize, y: usize) -> (f32, f32) {
        (self.u[x + y * self.width], self.v[x + y * self.width])
    }

//
// void FluidSolverJosStam::add_source(int N, float *&x, float *&s, float dt) {
//     int i, size = (N + 2) * (N + 2);
//     for ( i = 0; i < size; i++ ) {
//         x[i] += dt * s[i];
//     }
// }

    fn add_source(&mut self, N: usize, x: &mut Vec<f32>, s: &mut Vec<f32>, dt: f32) {
        for i in 0..(N + 2) * (N + 2) {
            x[i] += dt * s[i];
        }
    }

// void FluidSolverJosStam::set_bnd(int N, int b, float *&x) {
//     int i;
//     for ( i = 1; i <= N; i++ ) {
//         x[IX!(0, i)] = b == 1 ? -x[IX!(1, i)] : x[IX!(1, i)];
//         x[IX!(N + 1, i)] = b == 1 ? -x[IX!(N, i)] : x[IX!(N, i)];
//         x[IX!(i, 0)] = b == 2 ? -x[IX!(i, 1)] : x[IX!(i, 1)];
//         x[IX!(i, N + 1)] = b == 2 ? -x[IX!(i, N)] : x[IX!(i, N)];
//     }
//     x[IX!(0, 0)] = 0.5f * (x[IX!(1, 0)] + x[IX!(0, 1)]);
//     x[IX!(0, N + 1)] = 0.5f * (x[IX!(1, N + 1)] + x[IX!(0, N)]);
//     x[IX!(N + 1, 0)] = 0.5f * (x[IX!(N, 0)] + x[IX!(N + 1, 1)]);
//     x[IX!(N + 1, N + 1)] = 0.5f * (x[IX!(N, N + 1)] + x[IX!(N + 1, N)]);
// }

    fn set_bnd(&mut self, N: usize, b: usize, x: &mut Vec<f32>) {
        for i in 1..N {
            x[IX!(0, i)] = if b == 1 { -x[IX!(1, i)] } else { x[IX!(1, i)] };
            x[IX!(N + 1, i)] = if b == 1 { -x[IX!(N, i)] } else { x[IX!(N, i)] };
            x[IX!(i, 0)] = if b == 2 { -x[IX!(i, 1)] } else { x[IX!(i, 1)] };
            x[IX!(i, N + 1)] = if b == 2 { -x[IX!(i, N)] } else { x[IX!(i, N)] };
        }
        x[IX!(0, 0)] = 0.5f32 * (x[IX!(1, 0)] + x[IX!(0, 1)]);
        x[IX!(0, N + 1)] = 0.5f32 * (x[IX!(1, N + 1)] + x[IX!(0, N)]);
        x[IX!(N + 1, 0)] = 0.5f32 * (x[IX!(N, 0)] + x[IX!(N + 1, 1)]);
        x[IX!(N + 1, N + 1)] = 0.5f32 * (x[IX!(N, N + 1)] + x[IX!(N + 1, N)]);
    }

//
// void FluidSolverJosStam::lin_solve(int N, int b, float *&x, float *&x0, float a, float c) {
//     int i, j, k;
//
//     for ( k = 0; k < 2; k++ ) {
//         for ( i = 1; i <= N; i++ ) {
//             for ( j = 1; j <= N; j++ ) {
//                 x[IX!(i, j)] = (x0[IX!(i, j)] + a * (x[IX!(i - 1, j)] + x[IX!(i + 1, j)] + x[IX!(i, j - 1)] + x[IX!(i, j + 1)])) / c;
//             }
//         }
//         set_bnd(N, b, x);
//     }
// }

    fn lin_solve(&mut self, N: usize, b: usize, x: &mut Vec<f32>, x0: &mut Vec<f32>, a: f32, c: f32) {
        for k in 0..2 {
            for i in 1..N {
                for j in 1..N {
                    x[IX!(i, j)] = (x0[IX!(i, j)] + a * (x[IX!(i - 1, j)] + x[IX!(i + 1, j)] + x[IX!(i, j - 1)] + x[IX!(i, j + 1)])) / c;
                }
            }
            self.set_bnd(N, b, x);
        }
    }

//
// void FluidSolverJosStam::diffuse(int N, int b, float *&x, float *&x0, float viscosity, float dt) {
//     float a = dt * viscosity * N * N;
//     lin_solve(N, b, x, x0, a, 1 + 4 * a);
// }

    fn diffuse(&mut self, N: usize, b: usize, x: &mut Vec<f32>, x0: &mut Vec<f32>, viscosity: f32, dt: f32) {
        let a = dt * viscosity * N as f32 * N as f32;
        self.lin_solve(N, b, x, x0, a, 1.0f32 + 4.0f32 * a);
    }

//
// void FluidSolverJosStam::advect(int N, int b, float *&d, float *&d0, float *&u, float *&v, float dt) {
//     int i, j, i0, j0, i1, j1;
//     float x, y, s0, t0, s1, t1, dt0;
//
//     dt0 = dt * N;
//     for ( i = 1; i <= N; i++ ) {
//         for ( j = 1; j <= N; j++ ) {
//             x = i - dt0 * u[IX!(i, j)];
//             y = j - dt0 * v[IX!(i, j)];
//             if ( x < 0.5f ) {
//                 x = 0.5f;
//             }
//             if ( x > N + 0.5f ) {
//                 x = N + 0.5f;
//             }
//             i0 = (int) x;
//             i1 = i0 + 1;
//             if ( y < 0.5f ) {
//                 y = 0.5f;
//             }
//             if ( y > N + 0.5f ) {
//                 y = N + 0.5f;
//             }
//             j0 = (int) y;
//             j1 = j0 + 1;
//             s1 = x - i0;
//             s0 = 1 - s1;
//             t1 = y - j0;
//             t0 = 1 - t1;
//             d[IX!(i, j)] = s0 * (t0 * d0[IX!(i0, j0)] + t1 * d0[IX!(i0, j1)]) +
//                           s1 * (t0 * d0[IX!(i1, j0)] + t1 * d0[IX!(i1, j1)]);
//         }
//     }
//     set_bnd(N, b, d);
// }

    fn advect(&mut self, N: usize, b: usize, d: &mut Vec<f32>, d0: &mut Vec<f32>, u: &mut Vec<f32>, v: &mut Vec<f32>, dt: f32) {
        let dt0 = dt * N as f32;
        for i in 1..N {
            for j in 1..N {
                let mut x = i as f32 - dt0 * u[IX!(i, j)];
                let mut y = j as f32 - dt0 * v[IX!(i, j)];
                if x < 0.5f32 {
                    x = 0.5f32;
                }
                if x > N as f32 + 0.5f32 {
                    x = N as f32 + 0.5f32;
                }
                let i0 = x as usize;
                let i1 = i0 + 1;
                if y < 0.5f32 {
                    y = 0.5f32;
                }
                if y > N as f32 + 0.5f32 {
                    y = N as f32 + 0.5f32;
                }
                let j0 = y as usize;
                let j1 = j0 + 1;
                let s1 = x - i0 as f32;
                let s0 = 1.0f32 - s1;
                let t1 = y - j0 as f32;
                let t0 = 1.0f32 - t1;
                d[IX!(i, j)] = s0 * (t0 * d0[IX!(i0, j0)] + t1 * d0[IX!(i0, j1)]) +
                    s1 * (t0 * d0[IX!(i1, j0)] + t1 * d0[IX!(i1, j1)]);
            }
        }
        self.set_bnd(N, b, d);
    }

//
// void FluidSolverJosStam::project(int N, float *&u, float *&v, float *&p, float *&div) {
//     int i, j;
//
//     for ( i = 1; i <= N; i++ ) {
//         for ( j = 1; j <= N; j++ ) {
//             div[IX!(i, j)] = -0.5f * (u[IX!(i + 1, j)] - u[IX!(i - 1, j)] + v[IX!(i, j + 1)] - v[IX!(i, j - 1)]) / N;
//             p[IX!(i, j)] = 0;
//         }
//     }
//     set_bnd(N, 0, div);
//     set_bnd(N, 0, p);
//
//     lin_solve(N, 0, p, div, 1, 4);
//
//     for ( i = 1; i <= N; i++ ) {
//         for ( j = 1; j <= N; j++ ) {
//             u[IX!(i, j)] -= 0.5f * N * (p[IX!(i + 1, j)] - p[IX!(i - 1, j)]);
//             v[IX!(i, j)] -= 0.5f * N * (p[IX!(i, j + 1)] - p[IX!(i, j - 1)]);
//         }
//     }
//     set_bnd(N, 1, u);
//     set_bnd(N, 2, v);
// }

    fn project(&mut self, N: usize, u: &mut Vec<f32>, v: &mut Vec<f32>, p: &mut Vec<f32>, div: &mut Vec<f32>) {
        for i in 1..N {
            for j in 1..N {
                div[IX!(i, j)] = -0.5f32 * (u[IX!(i + 1, j)] - u[IX!(i - 1, j)] + v[IX!(i, j + 1)] - v[IX!(i, j - 1)]) / N as f32;
                p[IX!(i, j)] = 0.0f32;
            }
        }
        self.set_bnd(N, 0, div);
        self.set_bnd(N, 0, p);
        self.lin_solve(N, 0, p, div, 1.0f32, 4.0f32);
        for i in 1..N {
            for j in 1..N {
                u[IX!(i, j)] -= 0.5f32 * N as f32 * (p[IX!(i + 1, j)] - p[IX!(i - 1, j)]);
                v[IX!(i, j)] -= 0.5f32 * N as f32 * (p[IX!(i, j + 1)] - p[IX!(i, j - 1)]);
            }
        }
        self.set_bnd(N, 1, u);
        self.set_bnd(N, 2, v);
    }

//
// void FluidSolverJosStam::dens_step(int N, float *&x, float *&x0, float *&u, float *&v, float diff, float dt) {
//     add_source(N, x, x0, dt);
//     SWAP (x0, x);
//     diffuse(N, 0, x, x0, diff, dt);
//     SWAP (x0, x);
//     advect(N, 0, x, x0, u, v, dt);
// }
//

    fn dens_step(&mut self, N: usize, x: &mut Vec<f32>, x0: &mut Vec<f32>, u: &mut Vec<f32>, v: &mut Vec<f32>, diff: f32, dt: f32) {
        self.add_source(N, x, x0, dt);
        self.swap(x0, x);
        self.diffuse(N, 0, x, x0, diff, dt);
        self.swap(x0, x);
        self.advect(N, 0, x, x0, u, v, dt);
    }

// void FluidSolverJosStam::vel_step(int N, float *&u, float *&v, float *&u0, float *&v0, float viscosity, float dt) {
//     add_source(N, u, u0, dt);
//     add_source(N, v, v0, dt);
//     SWAP (u0, u);
//     diffuse(N, 1, u, u0, viscosity, dt);
//     SWAP (v0, v);
//     diffuse(N, 2, v, v0, viscosity, dt);
//     project(N, u, v, u0, v0);
//     SWAP (u0, u);
//     SWAP (v0, v);
//     advect(N, 1, u, u0, u0, v0, dt);
//     advect(N, 2, v, v0, u0, v0, dt);
//     project(N, u, v, u0, v0);
// }

    fn vel_step(&mut self, N: usize, u: &mut Vec<f32>, v: &mut Vec<f32>, u0: &mut Vec<f32>, v0: &mut Vec<f32>, viscosity: f32, dt: f32) {
        self.add_source(N, u, u0, dt);
        self.add_source(N, v, v0, dt);
        self.swap(u0, u);
        self.diffuse(N, 1, u, u0, viscosity, dt);
        self.swap(v0, v);
        self.diffuse(N, 2, v, v0, viscosity, dt);
        self.project(N, u, v, u0, v0);
        self.swap(u0, u);
        self.swap(v0, v);
        self.advect(N, 1, u, u0, u0, v0, dt);
        self.advect(N, 2, v, v0, u0, v0, dt);
        self.project(N, u, v, u0, v0);
    }

//
// void FluidSolverJosStam::printVelocity(){
//     std::cout << "u:" << std::endl;
//     for (int i = 1; i <= this->width; ++i) {
//         for (int j = 1; j <= this->height; ++j) {
//             std::cout << "(" << getVelocityU(i,j) << "," << getVelocityV(i,j) << ")";
//         }
//         std::cout << std::endl;
//     }
// }
//

    fn print_velocity(&self) {
        println!("u:");
        for i in 1..self.width {
            for j in 1..self.height {
                println!("({}, {})", self.get_velocity_u(i, j), self.get_velocity_v(i, j));
            }
            println!("");
        }
    }

    // void FluidSolverJosStam::fill(int direction){
//     int N = width;
//     for ( int y = 10; y < height-10; y++) {
//         if(std::cos(((float)y)*cos0)>cos1){
//             this->u0[IX!(3, y)] = advectAmount;
//         }
//     }
// }
    fn fill(&mut self, direction: usize) {
        let N = self.width;
        for y in 10..self.height - 10 {
            if (y as f32 * self.cos0).cos() > self.cos1 {
                self.u0[IX!(3, y)] = self.advect_amount;
            }
        }
    }

    //
// void FluidSolverJosStam::step(){
//     vel_step(width, u, v, u0, v0, viscosity, dt);
//     dens_step (width, dens, dens_prev, u, v, density, dt );
// }
    fn step(&mut self) {
        self.vel_step(self.width, &mut self.u, &mut self.v, &mut self.u0, &mut self.v0, self.viscosity, self.dt);
        self.dens_step(self.width, &mut self.dens, &mut self.dens_prev, &mut self.u, &mut self.v, self.density, self.dt);
    }

//
// float FluidSolverJosStam::getVelocityU(int x, int y){
//     int N = width;
//     return u[IX!(x,y)];
// }

    fn get_velocity_u(&self, x: usize, y: usize) -> f32 {
        let N = self.width;
        self.u[IX!(x, y)]
    }

    // float FluidSolverJosStam::getVelocityV(int x, int y){
//     int N = width;
//     return v[IX!(x,y)];
// }
    fn get_velocity_v(&self, x: usize, y: usize) -> f32 {
        let N = self.width;
        self.v[IX!(x, y)]
    }

//
// float FluidSolverJosStam::getPressure(int x, int y){
//     int N = width;
//     return dens[IX!(x,y)];
// }

    fn get_pressure(&self, x: usize, y: usize) -> f32 {
        let N = self.width;
        self.dens[IX!(x, y)]
    }

//
// FluidSolverJosStam::FluidSolverJosStam(int width, int height) {
//     this->width = width;
//     this->height = height;
//     int size = (width+2)*(height+2);
//     u			= (float *) malloc ( size*sizeof(float) );
//     v			= (float *) malloc ( size*sizeof(float) );
//     u0		    = (float *) malloc ( size*sizeof(float) );
//     v0		    = (float *) malloc ( size*sizeof(float) );
//     dens		= (float *) malloc ( size*sizeof(float) );
//     dens_prev	= (float *) malloc ( size*sizeof(float) );
//     int N = width;
// }

    fn new(width: usize, height: usize) -> Box<dyn Solver> {
        let size = (width + 2) * (height + 2);
        let u = vec![0.0; size];
        let v = vec![0.0; size];
        let u0 = vec![0.0; size];
        let v0 = vec![0.0; size];
        let dens = vec![0.0; size];
        let dens_prev = vec![0.0; size];
        let N = width;
        Solver {
            width: width,
            height: height,
            size: size,
            u: u,
            v: v,
            u0: u0,
            v0: v0,
            dens: dens,
            dens_prev: dens_prev,
            N: N,
            dt: 0.1,
            diff: 0.0,
            visc: 0.0,
            force: 5.0,
            source: 100.0,
            density: 0.0,
            viscosity: 0.0,
            cos0: 0.0,
            cos1: 0.0,
            advect_amount: 0.0,
        }
    }

//
// FluidSolverJosStam::~FluidSolverJosStam() {
//     if ( u ) free ( u );
//     if ( v ) free ( v );
//     if ( u0 ) free ( u0 );
//     if ( v0 ) free ( v0 );
//     if ( dens ) free ( dens );
//     if ( dens_prev ) free ( dens_prev );
// }

    fn drop(&mut self) {
        // if ( u ) free ( u );
        // if ( v ) free ( v );
        // if ( u0 ) free ( u0 );
        // if ( v0 ) free ( v0 );
        // if ( dens ) free ( dens );
        // if ( dens_prev ) free ( dens_prev );
    }
}