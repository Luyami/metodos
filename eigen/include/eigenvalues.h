#ifndef EIGENVALUES_H
#define EIGENVALUES_H

#include <mat.h>
#include <utility>

namespace EigenValues{
    static int MAX_ITERATIONS = 1000;

    pair<double, mat> compute_power_method(mat matrix, mat vec, double error);
    pair<double, mat> compute_inverse_power_method(mat matrix, mat vec, double error);
    pair<double, mat> compute_nearest_power_method(mat matrix, mat vec, double u, double error);
};

#endif