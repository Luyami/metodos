#include <transformations.h>

int main(){
    mat mat1(10,10);

    //Matrix to test the methods I implemented (: (Im only testing QR because QR uses householder and all the other algorithms implemented here. If one fails, QR must fail as well)
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

    std::pair<mat, mat> qr = Transformations::qrEigen(mat1);
    std::cout << "Eigen values matrix..." << '\n';
    qr.first.print();
    std::cout << '\n';

    std::cout << "Eigen vectors matrix..." << '\n';
    qr.second.print();
    std::cout << '\n';

    return 1;
}