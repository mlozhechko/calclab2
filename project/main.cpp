#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "system_solver.h"
#include "matrix_utils.h"

namespace ub = boost::numeric::ublas;

template <class T>
int lab2Main(int argc, char** argv) {
  bool isMatrixSpecified = cmdOptionExists(argv, argv + argc, "-matrix");
  bool isVectorSpecified = cmdOptionExists(argv, argv + argc, "-vector");
  if (!isMatrixSpecified || !isVectorSpecified) {
    std::cerr << "source data is not specified" << std::endl;
    return -1;
  }

  if (cmdOptionExists(argv, argv + argc, "-debug")) {
    std::cout << "debug mode enabled" << std::endl;
    log::setEnabled(true);
    log::setPrecision(5);
  }

  if (!cmdOptionExists(argv, argv + argc, "-method")) {
    std::cerr << "solver method is not specified" << std::endl;
    return -2;
  }

  std::ifstream matrixStream(getCmdOption(argv, argv + argc, "-matrix"));
  std::ifstream vectorStream(getCmdOption(argv, argv + argc, "-vector"));

  ub::matrix<T> A, B;
  initMatrix(A, matrixStream);
  initVector(B, vectorStream);

  std::cout << "source A: " << A << std::endl;
  std::cout << "source B: " << B << std::endl;

  ub::matrix<T> X;
  std::string method(getCmdOption(argv, argv + argc, "-method"));

  /*
   * calculation
   */
  T eps = 0.001;
  std::function<T(const ub::matrix<T>&)> norm = normCubic<T>;

  if (method == "fpi") {
    std::cout << "using fixed point iteration method" << std::endl;
    /*
     * fpi works correctly with 6 test with 0.25 tau and Cubic norm
     * with this conditions norm(C) <= 1 (0.95)
     */
    T tau = 0.05;
    if (fixedPointIteration(A, B, X, tau, norm, eps) < 0) {
      return -5;
    }
  } else if (method == "jacobi") {
    std::cout << "using jacobi method" << std::endl;
    if (jacobiIteration(A, B, X, norm, eps) < 0) {
      return -1;
    }
  } else if (method == "seidel") {
    std::cout << "using seidel method" << std::endl;
    if (zeidelIteration(A, B, X, norm, eps) < 0) {
      return -1;
    }
  } else if (method == "relax3d") {
    std::cout << "using relaxation method (3-diagonal matrices case)" << std::endl;
    std::cout << "warning(!) source matrices will be redefined" << std::endl;

    const T w = 1;
    const ssize_t N = 213;
    A = ub::zero_matrix<T>(N, 3);
    B = ub::zero_matrix<T>(N, 1);

    A(0, 1) = 4;
    A(0, 2) = 1;

    for (ssize_t i = 0; i < N - 1; ++i) {
      A(i, 0) = 1;
      A(i, 1) = 4;
      A(i, 2) = 1;
    }
    A(N - 1, 0) = 1;
    A(N - 1, 1) = 4;

    B(0, 0) = 6;
    for (ssize_t i = 1; i < N - 1; ++i) {
      B(i, 0) = 10 - 2 * ((i + 1) % 2);
    }
    B(N - 1, 0) = 9 - 3 * (N % 2);

    if (diag3RelaxaionIteration(A, B, X, norm, eps, w) < 0) {
      return -1;
    }
  } else {
      std::cerr << "solver method cannot be parsed" << std::endl;
      return -3;
  }

  std::cout << "result is X = " << X << std::endl;

  return 0;
}


int main(int argc, char** argv) {
  if (!cmdOptionExists(argv, argv + argc, "-precision")) {
    std::cerr << "precision is not specified" << std::endl;
    return -1;
  }

  std::string precision = getCmdOption(argv, argv + argc, "-precision");
  if (precision == "double") {
    return lab2Main<double>(argc, argv);
  } else if (precision == "float") {
    return lab2Main<float>(argc, argv);
  }

  std::cerr << "precision cannot be parsed correctly" << std::endl;
  return -2;
}
