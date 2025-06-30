#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <iostream>

//vec2 when z = 0
struct vec3{
    double x;
    double y;
    double z;

    double norm();
    void print();

    static vec3 add(const vec3& v1, const vec3& v2);
    static vec3 sub(const vec3& v1, const vec3& v2);
    static vec3 mult(const vec3& v, double n);
    static vec3 div(const vec3& v, double n);
    static double dot(const vec3& v1, const vec3& v2);
};

#endif