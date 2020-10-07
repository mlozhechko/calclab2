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
  if (method == "fpi") {
    std::cout << "using fixed point iteration method" << std::endl;

    T tau = 1;
    T eps = std::numeric_limits<T>::epsilon();
    if (fixedPointIteration(A, B, X, tau, 0, eps) < 0) {
      return -5;
    }
  } else if (method == "jacobi") {
    std::cout << "using jacobi method" << std::endl;

    T eps = 0.1;
    std::function<T(const ub::matrix<T>&)> norm = normOcta<T>;
    if (jacobiIteration(A, B, X, norm, eps) < 0) {
      return -1;
    }
  } else {
    std::cerr << "solver method cannot be parsed" << std::endl;
    return -3;
  }

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
