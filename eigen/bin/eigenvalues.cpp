#include <eigenvalues.h>

pair<double, mat> EigenValues::compute_power_method(mat matrix, mat vec, double error){
    if (matrix.n != matrix.m) throw invalid_argument("Must be a square matrix!");
    if (matrix.n != vec.n)    throw invalid_argument("Vector rows must match square matrix's order in power method!");
    if (vec.m > 1)            throw invalid_argument("Variable vec is not a column vector!");

    int iterations = 0;

    double eig_val_old = 0;
    double eig_val_new = 100;

    mat old_vec = vec;
    mat new_vec = vec;

    while ( abs((eig_val_new - eig_val_old)/eig_val_new) > error ){
        old_vec = new_vec;
        new_vec = mat::mult(matrix, old_vec);
        new_vec = mat::normalize(new_vec);

        eig_val_old = eig_val_new;
        eig_val_new = (mat::transpose(new_vec) * matrix * new_vec).matrix[0][0]/
                      (mat::transpose(new_vec) * new_vec)         .matrix[0][0];

        ++iterations;

        if (iterations >= EigenValues::MAX_ITERATIONS) throw runtime_error("Convergence issues in power method!");
    }

    return pair(eig_val_new, old_vec);
}

pair<double, mat> EigenValues::compute_inverse_power_method(mat matrix, mat vec, double error){
    if (matrix.n != matrix.m) throw invalid_argument("Must be a square matrix!");

    mat inverse = mat::inverse(matrix);
    
    pair<double, mat> values = EigenValues::compute_power_method(inverse, vec, error);

    return pair(1.0f/values.first, values.second);
}

pair<double, mat> EigenValues::compute_nearest_power_method(mat matrix, mat vec, double u, double error){
    if (matrix.n != matrix.m) throw invalid_argument("Must be a square matrix!");

    mat nearestMat = matrix - mat::identity(matrix.n) * u;

    pair<double, mat> values = EigenValues::compute_inverse_power_method(nearestMat, vec, error);

    return pair(values.first + u, values.second);
}