/**
 * Just a simple copy'n'paste template class to quickly add tests
 */

#include "Base.hh"

#include "../src/solver.hpp"

static const int TEST_ITERATIONS = 3;
using namespace std;

/**
 * Test usage, shape and correct indexing of solver and results produced by eigen library.
 */
void test_eigen_linear_solver_shape() {
    int width = 3;
    int height = 3;
    int square_size = width * height;

    SparseMatrix<double> A(square_size, square_size);

    VectorXd rhs(square_size);
    VectorXd p(square_size);

    /* Matrix A is sparse. If we rawly know the number of non zero coefficients per row, we can tell the matrix
     * beforehand to pre-allocate the according number of columns. For big matrices this speeds up solving
     * significantly!  */
    A.reserve(VectorXd::Constant( square_size,6));

    int i=0;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            i = x + y*width;
            A.coeffRef(i,i) = i+1;
            rhs.coeffRef(i) = square_size-i;
        }
    }

    std::cout << "A = " << std::endl;

    for (int x = 0; x < square_size; x++) {
        for (int y = 0; y < square_size; y++) {
            std::cout << A.coeff(x,y) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << "rhs = " << std::endl;

    int x=0,y=0;
    for (int idx=0; idx < square_size; idx++) {
        x = idx % width;
        y = idx / width;
        std::cout << rhs(idx) << " // "<<  x << "," << y << std::endl;
    }

    Eigen::ConjugateGradient<SparseMatrix<double>> solver(A);

    p = solver.solve(rhs);

    std::cout << std::endl;

    std::cout << "p = " << std::endl;

    for (int idx=0; idx < square_size; idx++) {
        std::cout << p(idx) << std::endl;
    }

}

int main(int argc, const char *argv[]) {
    test_eigen_linear_solver_shape();
    return 0;
}