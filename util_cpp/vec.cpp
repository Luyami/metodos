#include <vec.h>

double vec3::norm(){
    return std::sqrt(
        this->x*this->x + this->y*this->y + this->z*this->z
    );
}

void vec3::print(){
    std::cout << "(" << x << ", " << y << ", " << z << ")" << '\n';
}

vec3 vec3::add(const vec3& v1, const vec3& v2) {
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

vec3 vec3::sub(const vec3& v1, const vec3& v2) {
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

vec3 vec3::mult(const vec3& v, double n) {
    return { v.x * n, v.y * n, v.z * n };
}

vec3 vec3::div(const vec3& v, double n) {
    return { v.x / n, v.y / n, v.z / n };
}

double vec3::dot(const vec3& v1, const vec3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}