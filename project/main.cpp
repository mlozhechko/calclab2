#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "system_solver.h"
#include "matrix_utils.h"

namespace ub = boost::numeric::ublas;

template <class T>
int lab1Main(int argc, char** argv) {
  /*
   * initialization
   */
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
  if (method == "qr") {
    std::cout << "using qr method" << std::endl;
    if (QRSolve(A, B, X) < 0) {
      return -5;
    }
  } else if (method == "gauss") {
    std::cout << "using gauss method" << std::endl;
    if (gaussSolve(A, B, X) < 0) {
      return -4;
    }
  } else {
    std::cerr << "solver method cannot be parsed" << std::endl;
    return -3;
  }
  std::cout << "X: " << X << std::endl;

  ub::matrix<T> Bm;
  matrixMult(A, X, Bm);

  ub::matrix<T> deltaB = Bm - B;

  log::debug() << "calculated B: " << Bm << "\n";

  std::cout << "residual Octa norm: " << normOcta(deltaB) << std::endl;
  std::cout << "residual Cubic norm: " << normCubic(deltaB) << std::endl;


  ub::matrix<T> AI;
  invertMatrix(A, AI);
  std::cout << "cubic cond A = " << normCubic(AI) * normCubic(A) << std::endl;
  std::cout << "octa cond A = " << normOcta(AI) * normOcta(A) << std::endl;

  ub::matrix<T> AM;
  matrixMult(A, AI, AM);
  std::cout << "A * AInv = " << AM << std::endl;

  if (!cmdOptionExists(argv, argv + argc, "-pert")) {
    return 0;
  }

  std::cout << "perturbations module. cond estimation" << std::endl;
  /*
   * system perturbations
   */
  size_t counter = 0;
  for (auto pert : {-0.1, -0.05, 0.05, 0.1}) {
    for (ssize_t i = 0; i < B.size1(); ++i) {
      ub::matrix<T> BP = B;
      ub::matrix<T> XP;
      BP(i, 0) += pert;

      if (method == "qr") {
        QRSolve(A, BP, XP);
      } else if (method == "gauss") {
        gaussSolve(A, BP, XP);
      } else {
        std::cerr << "method cannot be specified" << std::endl;
        return -4;
      }

      ub::matrix<T> deltaX = X - XP;
      ub::matrix<T> deltaBP = B - BP;

      std::cout << "[pert iteration " << ++counter << "]" << std::endl;
      std::cout << "B perturbation row = " << i << " with value = " << pert << std::endl;
      std::cout << "original solution: " << X << std::endl;
      std::cout << "perturbed solution: " << XP << std::endl;
      log::debug() << "solutions delta: " << deltaX << "\n";

      T bx = 0, bb = 0;
      bx = normOcta(deltaX) / normOcta(X);
      bb = normOcta(deltaBP) / normOcta(B);
      std::cout << "octa condA >= " << bx / bb << std::endl;

      bx = normCubic(deltaX) / normCubic(X);
      bb = normCubic(deltaBP) / normCubic(B);
      std::cout << "cubic condA >= " << bx / bb << std::endl;
    }
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
    return lab1Main<double>(argc, argv);
  } else if (precision == "float") {
    return lab1Main<float>(argc, argv);
  }

  std::cerr << "precision cannot be parsed correctly" << std::endl;
  return -2;
}
