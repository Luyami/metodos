#include <mat.h>
#include <vec.h>
#include <eigenvalues.h>
#include <iostream>

static double ERROR = 0.00000000000001;

int main(){
    mat mat1 = mat(10,10);
    mat vec1 = mat(10,1);

    //Matrix to test the methods I implemented (:
    mat1.build({
    0.62, 0.49, 0.56, 0.42, 0.45, 0.31, 0.59, 0.67, 0.52, 0.60,
    0.49, 0.97, 0.65, 0.71, 0.58, 0.62, 0.69, 0.54, 0.50, 0.64,
    0.56, 0.65, 0.71, 0.68, 0.60, 0.63, 0.55, 0.53, 0.59, 0.57,
    0.42, 0.71, 0.68, 0.89, 0.67, 0.64, 0.61, 0.66, 0.62, 0.63,
    0.45, 0.58, 0.60, 0.67, 0.79, 0.55, 0.58, 0.60, 0.57, 0.52,
    0.31, 0.62, 0.63, 0.64, 0.55, 0.83, 0.65, 0.59, 0.60, 0.61,
    0.59, 0.69, 0.55, 0.61, 0.58, 0.65, 0.91, 0.66, 0.64, 0.67,
    0.67, 0.54, 0.53, 0.66, 0.60, 0.59, 0.66, 0.87, 0.62, 0.55,
    0.52, 0.50, 0.59, 0.62, 0.57, 0.60, 0.64, 0.62, 0.73, 0.63,
    0.60, 0.64, 0.57, 0.63, 0.52, 0.61, 0.67, 0.55, 0.63, 0.81
    });

    vec1.build({16.3,73.2,20.1,12.05,5.03,10.2,34.2,100.2,4.2,30.1});

    std::cout << '\n';

    std::cout << "Power method..." << '\n';
    pair<double, mat> values1 = EigenValues::compute_power_method(mat1, vec1, ERROR);
    std::cout << setprecision(10) << "Eigen Value: " << values1.first << '\n';
    values1.second.print();
    std::cout << "\n";

    std::cout << "Inverse power method..." << '\n';
    pair<double, mat> values2 = EigenValues::compute_inverse_power_method(mat1, vec1, ERROR);
    std::cout << setprecision(10) << "Eigen Value: " << values2.first << '\n';
    values2.second.print();
    std::cout << "\n";

    std::cout << "Nearest power method..." << '\n';
    pair<double, mat> values3 = EigenValues::compute_nearest_power_method(mat1, vec1, 0.4, ERROR);
    std::cout << setprecision(10) << "Eigen Value: " << values3.first << '\n';
    values3.second.print();
    std::cout << "\n";

    std::cout << '\n';

    return 0;
}