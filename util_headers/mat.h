#ifndef MAT_H
#define MAT_H

#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>

using namespace std;

struct mat;

//Square matrix of size nxm
struct mat{
    mat(uint32_t size); //For square mat where n = m
    mat(uint32_t n, uint32_t m);
    mat(const mat& other);

    uint32_t n;
    uint32_t m;
    double** matrix;

    void build(vector<double> values);
    void print() const;

    double norm() const; //2-Norm (euclidian)

    static mat add(const mat& m1, const mat& m2);
    static mat sub(const mat& m1, const mat& m2);
    static mat mult(const mat& m1, const mat& m2);
    static mat mult(const mat& m, double scalar);
    static mat div(const mat& m, double scalar);
    static mat normalize(const mat& m);
    static mat transpose(const mat& m);
    static mat inverse(const mat& m);
    static mat identity(uint32_t size);
    static double det(const mat& m);
    static bool isSymmetric(const mat& m);

    mat operator+(const mat& other) const;
    mat operator-(const mat& other) const;
    mat operator*(const mat& other) const;    
    mat operator*(double scalar) const;       
    mat operator/(double scalar) const;        
    mat& operator=(const mat& other);

    ~mat();
};

#endif