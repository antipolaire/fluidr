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
#include <Eigen/IterativeLinearSolvers>
// #include <Eigen/CholmodSupport>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>

#include "../src/bnas.h"
#include "../src/field2d.hpp"
#include "../src/MACGrid.hpp"

#define PRINT_DATA_INFOS

using namespace Eigen;
using namespace std;

void printLGS(const SparseMatrix<double> &A, const VectorXd &b, const VectorXd &x);

void
test_sparse_solver_performance(const std::function<void(SparseMatrix<double> &, VectorXd &, VectorXd &)> &fnSolver) {
    /* Example how to solve Ax=b using conjugated gradient algorithm with incomplete cholesky preconditioner using eigen library */

    int n = 10;

    // Square (Coefficient-) Matrix
    SparseMatrix<double> A = SparseMatrix<double>(n, n);

    // Right hand side known column vector b
    VectorXd b(n);

    /* Unknown column vector x*/
    VectorXd x(n);

    VectorXd x_test(n);

    x_test(0) = 0.133966;
    x_test(1) = 0.535863;
    x_test(2) = 1.00949;
    x_test(3) = 1.50208;
    x_test(4) = 1.99884;
    x_test(5) = 2.49328;
    x_test(6) = 2.97428;
    x_test(7) = 3.40382;
    x_test(8) = 3.64102;
    x_test(9) = 3.16025;

    /**/
    for (int i = 0; i < n; i++) {
        b(i) = i;
        A.coeffRef(i, i) = 4;
        if (i < n - 1) {
            A.coeffRef(i + 1, i) = -1;
            A.coeffRef(i, i + 1) = -1;
        }
    }

    /* Print what we have initially */
    // printLGS(A, b, x);

    std::clock_t c_start = std::clock();
    /* Read http://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html for other solvers */

    /* Solve it */
    // x = solver.solve(b);

    fnSolver(A, x, b);

    std::clock_t c_end = std::clock();

    long time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "cpu time: " << time_elapsed_ms << "ms" <<
              std::endl;

    /* Print what the solver gave us */
    // printLGS(A, b, x);

    for (int i = 0; i < n; i++) {
        if (x(i) != x_test(i)) {
            cout << "unexpected solution:" << std::endl;

            cout << "expected x=" << std::endl << x_test;
            cout << std::endl;
            cout << "actual x=" << std::endl << x;
            cout << std::endl;
            break;
        }
    }
}

void test_all_solvers_performance() {

    std::cout << "Direct sparse LDLT Cholesky factorizations without square root:" << std::endl;

    test_sparse_solver_performance([](SparseMatrix<double> &A, VectorXd &x, VectorXd &b) {
        Eigen::SimplicialLDLT<SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver(A);
        x = solver.solve(b);
    });

    std::cout << "Conjugate gradient solver using incomplete cholesky preconditioning" << std::endl;

    test_sparse_solver_performance([](SparseMatrix<double> &A, VectorXd &x, VectorXd &b) {
        Eigen::ConjugateGradient<SparseMatrix<double>,
                Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> solver(A);
        x = solver.solve(b);
    });

    std::cout << "Direct sparse LLT Cholesky factorization" << std::endl;

    test_sparse_solver_performance([](SparseMatrix<double> &A, VectorXd &x, VectorXd &b) {
        Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver(A);
        x = solver.solve(b);
    });

    std::cout << "Least square Conjugate gradient solver" << std::endl;

    test_sparse_solver_performance([](SparseMatrix<double> &A, VectorXd &x, VectorXd &b) {
        Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>> solver(A);
        x = solver.solve(b);
    });


    std::cout << "Bi conjugate gradient stabilized solver" << std::endl;

    test_sparse_solver_performance([](SparseMatrix<double> &A, VectorXd &x, VectorXd &b) {
        Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver(A);
        x = solver.solve(b);
    });

    /*
    std::cout << "Cholmod simplicial direct cholesky (LLT) solver " << std::endl;

    test_sparse_solver_performance([](SparseMatrix<double> &A, VectorXd &x, VectorXd &b) {
        Eigen::CholmodSimplicialLLT<SparseMatrix<double>> solver(A);
        x = solver.solve(b);
    });

    std::cout << "Cholmod supernodal cholesky (LLT) solver " << std::endl;

    test_sparse_solver_performance([](SparseMatrix<double> &A, VectorXd &x, VectorXd &b) {
        Eigen::CholmodSupernodalLLT<SparseMatrix<double>> solver(A);
        x = solver.solve(b);
    });
    */

}

int main(int argc, const char *argv[]) {
    test_all_solvers_performance();
    return 0;
}

void printLGS(const SparseMatrix<double> &A, const VectorXd &b, const VectorXd &x) {
#ifdef PRINT_DATA_INFOS
    cout << "A: " << endl << static_cast<const SparseMatrixBase<SparseMatrix<double>> &>(A) << endl;
    cout << "b: " << endl << b << endl;
    cout << "x: " << endl << x << endl;
#endif
}