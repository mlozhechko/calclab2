#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <numeric>
#include "utils.h"

namespace ub = boost::numeric::ublas;

/*
 * synopsis:
 */
template<class T>
int gaussSolve(const ub::matrix<T>& sourceMatrix, const ub::matrix<T>& sourceVector,
               ub::matrix<T>& solution);

template <class T>
int initGivensCoefficients(const ub::matrix<T>& sourceMatrix, ssize_t row1, ssize_t row2,
                                  T& c, T& s);

template <class T>
int matrixTranspose(ub::matrix<T>& A);

template <class T>
int matrixFastRotate(ub::matrix<T>& A, ssize_t row1, ssize_t row2, T c, T s);

template<class U>
int QRDecomposition(const ub::matrix<U>& sourceMatrix, ub::matrix<U>& Q, ub::matrix<U>& R);

template<class T>
int QRSolve(const ub::matrix<T>& sourceMatrix, const ub::matrix<T>& sourceVector,
            ub::matrix<T>& solution);

template <class T>
int matrixMult(const ub::matrix<T>& sourceMatrixA, const ub::matrix<T>& sourceMatrixB,
               ub::matrix<T>& result);

template <class T>
int invertMatrix(const ub::matrix<T>& A, ub::matrix<T>& result);
/*
 * implementation:
 */
template<class T>
int gaussSolve(const ub::matrix<T>& sourceMatrix, const ub::matrix<T>& sourceVector,
               ub::matrix<T>& solution) {
  ub::matrix<T> A = sourceMatrix;
  ub::matrix<T> B = sourceVector;
  auto& X = solution;

  /*
   * ublas::matrix using linear memory model
   * therefore we cannot simply swap rows
   * and have to use transpositions matrix
   */
  ssize_t n = A.size1();
  std::vector<ssize_t> transposeVec(n);
  std::iota(transposeVec.begin(), transposeVec.end(), 0);
  auto& TR = transposeVec;

  for (ssize_t k = 0; k < n - 1; ++k) {
    T columnAbsMax = std::abs(A(TR[k], k));
    ssize_t columnAbsMaxIndex = k;
    for (ssize_t i = k + 1; i < n; ++i) {
      T absValue = std::abs(A(TR[i], k));
      if (absValue > columnAbsMax) {
        columnAbsMax = absValue;
        columnAbsMaxIndex = i;
      }
    }

    if (k != columnAbsMaxIndex) {
      std::swap(transposeVec[k], transposeVec[columnAbsMaxIndex]);
    }

    if (std::abs(A(TR[k], k)) < std::numeric_limits<T>::epsilon()) {
      std::cerr << "matrix det(A) == 0, system has infinite set of solutions" << std::endl;
      return -1;
    }

    for (ssize_t i = k + 1; i < n; ++i) {
      T c = A(TR[i], k) / A(TR[k], k);

      if (log::isDebug()) {
          A(TR[i], k) = 0;
      }

      for (ssize_t j = k + 1; j < n; ++j) {
        A(TR[i], j) -= c * A(TR[k], j);
      }
      B(TR[i], 0) -= c * B(TR[k], 0);

    }

    auto& logger = log::debug();
    logger << "cycle: " << k << "; TR: ";
    for (const auto& it: transposeVec) {
      logger << it << " ";
    }
    logger << "\n";
    log::debug() << "after gauss A: " << A << "\n";
    log::debug() << "after gauss B: " << B << "\n";
  }

  if (std::abs(A(TR[n - 1], n - 1)) < std::numeric_limits<T>::epsilon()) {
    std::cerr << "matrix det(A) == 0, system has infinite set of solutions" << std::endl;
    return -1;
  }

  auto& logger = log::debug();
  logger << "after gauss TR: ";
  for (const auto& it: transposeVec) {
    logger << it << " ";
  }
  logger << "\n";
  log::debug() << "after gauss A: " << A << "\n";
  log::debug() << "after gauss B: " << B << "\n";

  X.resize(n, 1);
  for (ssize_t i = n - 1; i >= 0; --i) {
    T ASum = 0;
    for (ssize_t j = n - 1; j > i; --j) {
      ASum += A(TR[i], j) * X(j, 0);
    }

    if (std::abs(A(TR[i], i)) < std::numeric_limits<T>::epsilon()) {
      std::cerr << "gauss calculation error" << std::endl;
    }

    X(i, 0) = (B(TR[i], 0) - ASum) / A(TR[i], i);
  }

  log::debug() << "result: " << X << "\n";

  return 0;
}

template<class T>
int initGivensCoefficients(const ub::matrix<T>& sourceMatrix, ssize_t row1, ssize_t row2,
                                  T& c, T& s) {
  const ub::matrix<T>& A = sourceMatrix;

  T denom = A(row1, row1) * A(row1, row1) + A(row2, row1) * A(row2, row1);
  denom = std::sqrt(denom);

  if (denom < std::numeric_limits<T>::epsilon()) {
    return -1;
  }

  c = A(row1, row1) / denom;
  s = A(row2, row1) / denom;

  return 0;
}


template<class T>
int matrixTranspose(ub::matrix<T>& A) {
  if (A.size1() != A.size2()) {
    std::cerr << "transpose of non-square matrices is not supported" << std::endl;
    return -1;
  }

  const ssize_t height = A.size1();
  const ssize_t width = A.size2();
  for (ssize_t i = 0; i < height; ++i) {
    for (ssize_t j = i + 1; j < width; ++j) {
      std::swap(A(i, j), A(j, i));
    }
  }
  return 0;
}

template<class T>
int matrixFastRotate(ub::matrix<T>& A, ssize_t row1, ssize_t row2, T c, T s) {
  ssize_t width = A.size2();
  for (size_t i = 0; i < width; ++i) {
    T tmp = A(row1, i) * c + A(row2, i) * s;
    A(row2, i) = A(row1, i) * (-s) + A(row2, i) * c;
    A(row1, i) = tmp;
  }

  return 0;
}

template<class U>
int QRDecomposition(const ub::matrix<U>& sourceMatrix, ub::matrix<U>& Q, ub::matrix<U>& R) {
  const ssize_t n = sourceMatrix.size1();

  ub::matrix<U> T = ub::identity_matrix(n, n);
  R = sourceMatrix;

  for (ssize_t i = 0; i < n - 1; ++i) {
    bool isDegenerateMatrix = true;

    for (ssize_t j = i + 1; j < n; ++j) {
      U c = 0, s = 0;
      log::debug() << "QR i: " << i << " j: " << j << "\n";
      log::debug() << "Q matrix " << T << "\n";
      log::debug() << "R matrix " << R << "\n";

      if (initGivensCoefficients(R, i, j, c, s) >= 0) {
        isDegenerateMatrix = false;

        log::debug() << "performing rotation c: " << c << " s: " << s << "\n";
        matrixFastRotate(T, i, j, c, s);
        matrixFastRotate(R, i, j, c, s);
      }
    }

    if (isDegenerateMatrix) {
      std::cerr << "matrix det == 0, QR decomposition cannot be completed" << std::endl;
      return -1;
    }
  }

  if (std::abs(R(n - 1, n - 1)) < std::numeric_limits<U>::epsilon()) {
    std::cerr << "matrix det == 0, QR decomposition cannot be completed" << std::endl;
    return -1;
  }

  log::debug() << "Q matrix " << T << "\n";
  log::debug() << "R matrix " << R << "\n";

  if (log::isDebug()) {
    auto QT = T;
    matrixTranspose(QT);
    ub::matrix<U> I;
    matrixMult(T, QT, I);
    log::debug() << "QT Q * Q^T = " << I << "\n";
  }

  matrixTranspose(T);
  Q = std::move(T);
  return 0;
}

template <class T>
int QRSolve(const ub::matrix<T>& sourceMatrix, const ub::matrix<T>& sourceVector,
            ub::matrix<T>& solution) {
  ub::matrix<T> Q, R;
  if (QRDecomposition(sourceMatrix, Q, R) < 0) {
    std::cerr << "system cannot be solved with QR decomposition method" << std::endl;
    return -1;
  }

  log::debug() << "result of QR decomposition" << "\n";
  log::debug() << "Q: " << Q << "\n";
  log::debug() << "R: " << R << "\n";

  ub::matrix<T> Bs;
  matrixTranspose(Q);
  matrixMult(Q, sourceVector, Bs);

  ssize_t n = sourceVector.size1();
  ub::matrix<T>& X = solution;
  X.resize(n, 1);
  for (ssize_t i = n - 1; i >= 0; --i) {
    T ASum = 0;
    for (ssize_t j = n - 1; j > i; --j) {
      ASum += R(i, j) * X(j, 0);
    }

    if (std::abs(R(i, i)) < std::numeric_limits<T>::epsilon()) {
      std::cerr << "QR calculation error" << std::endl;
    }

    X(i, 0) = (Bs(i, 0) - ASum) / R(i, i);
  }

  return 0;
}

template <class T>
int matrixMult(const ub::matrix<T>& sourceMatrixA, const ub::matrix<T>& sourceMatrixB,
               ub::matrix<T>& result) {
  const ub::matrix<T>& A = sourceMatrixA;
  const ub::matrix<T>& B = sourceMatrixB;
  ub::matrix<T>& X = result;

  if (A.size2() != B.size1()) {
    log::debug() << "matrices cannot be multiplied" << "\n";
    return -1;
  }
  const ssize_t len = A.size2();
  const ssize_t height = A.size1();
  const ssize_t width = B.size2();

  X.resize(height, width);
  for (ssize_t i = 0; i < height; ++i) {
    for (ssize_t j = 0; j < width; ++j) {
      T val = 0;
      for (ssize_t k = 0; k < len; ++k) {
        val += A(i, k) * B(k, j);
      }
      X(i, j) = val;
    }
  }

  return 0;
}

template<class T>
int invertMatrix(const ub::matrix<T>& A, ub::matrix<T>& result) {
  if (A.size1() != A.size2()) {
    std::cerr << "non square matrices invert is not supported" << std::endl;
  }
  ssize_t n = A.size1();

  ub::matrix<T>& AI = result;
  AI = ub::zero_matrix<T>(n, n);

  ub::matrix<T> Q, R;
  QRDecomposition(A, Q, R);

  matrixTranspose(Q);

  for (ssize_t k = 0; k < n; ++k) {
    for (ssize_t i = n - 1; i >= 0; --i) {
      T ASum = 0;
      for (ssize_t j = n - 1; j > i; --j) {
        ASum += R(i, j) * AI(j, k);
      }

      if (std::abs(R(i, i)) < std::numeric_limits<T>::epsilon()) {
        std::cerr << "QR calculation error" << std::endl;
      }

      AI(i, k) = (Q(i, k) - ASum) / R(i, i);
    }
  }

  if (log::isDebug()) {
    ub::matrix<T> multRes;
    matrixMult(A, AI, multRes);

    log::debug() << "invert of matrix result: " << multRes << "\n";
  }

  return 0;
}