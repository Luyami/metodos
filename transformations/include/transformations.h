#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <mat.h>
#include <utility>

namespace Transformations{
    static int MAX_ITERATIONS = 200;

    mat getHouseholderForSymmetricMat(const mat& A, int i);
    mat tridiagonalize(const mat& A);
    pair<mat, mat> qrDecompose(const mat& A);
    pair<mat, mat> qrEigen(const mat& A);
}

#endif