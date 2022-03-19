//
// Created by indevn on 22.3.18.
//

#ifndef SMALLPT_MATH3D_HPP
#define SMALLPT_MATH3D_HPP

#include <cmath>

class Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z;                  // position, also color (r,g,b)
    Vec (double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+ (const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator- (const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator* (double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult (const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &normalize () { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot (const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vec operator% (Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};

class Vector3:Vec{};

class Point3:Vector3{};
class Ray:Vector3{};

//class math3d {
//
//};


#endif //SMALLPT_MATH3D_HPP
