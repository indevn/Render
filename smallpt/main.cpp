// http://www.kevinbeason.com/smallpt
#include "erand48.h"// pseudo-random number generator, not in <stdlib.h>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -origin smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2

#define M_PI 3.141592653589793238462643


struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z;                  // position, also color (r,g,b)
    Vec (double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+ (const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator- (const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator* (double b) const { return Vec(x * b, y * b, z * b); }
    Vec operator* (const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec mult (const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &normalize () { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot (const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vec operator% (Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};

struct Vector3 : Vec {};

struct Point3 : Vector3 {};


struct Ray {
    Vec origin, direction;
    Ray (Vec o_, Vec d_) : origin(o_), direction(d_) {}
};

// 材质类型，diffuse，specular，菲涅尔镜面散射fresnel specular scattering
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Sphere {
    double radius;       // radius
    Vec possion; // center position of Sphere
    Vec emission, color;      // 表面自发光亮度 emission, color
    Refl_t refl;      // 材质类型reflection type (DIFFuse, SPECular, REFRactive)
    Sphere (double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
            radius(rad_), possion(p_), emission(e_), color(c_), refl(refl_) {}

    // 光线与球体求交：本质是解一个一元二次方程
    double intersect (const Ray &r) const { // returns distance, 0 if nohit
        Vec op = possion - r.origin;
        // Solve t^2*direction.direction + 2*t*(origin-possion).direction + (origin-possion).(origin-possion)-R^2 = 0

        double eps = 1e-4;
        double b = op.dot(r.direction);
        double det = b * b - op.dot(op) + radius * radius;


        if (det < 0) // det<0 means no hit
            return 0;
        else
            det = sqrt(det);

        double t;

        // cpp17以前的写法
//        if ((t = b - det) > eps) // t>0 means hit
//            return t;
//        else if ((t = b + det) > eps) // t<0 means no hit
//            return t;
//        else
//            return 0;

        // cpp17特性：带初始化器的if语句
        if (t = b - det;t > eps) // t>0 means hit
            return t;
        else if (t = b + det;t > eps) // t<0 means no hit
            return t;
        else
            return 0;

        // 不易维护的写法
//        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

// 场景内容维护
Sphere spheres[] = {//Scene: radius, position, emission, color, material
        Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),//Left
        Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF),//Rght
        Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),//Back
        Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),//Frnt
        Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Botm
        Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Top
        Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),//Mirr
        Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),//Glas
        Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF) //Lite
};

inline double clamp (double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt (double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); } // Gamma Encoding 将线性空间计算出来的颜色保存到sRGB空间


inline bool intersect (const Ray &r,
                       double &minDistance, // distance to intersection
                       int &id) {

    double n = sizeof(spheres) / sizeof(Sphere);
    double inf = minDistance = 1e20;

    double distance;

    // 逗号运算符：二元，先求解表达式 1，再求解表达式，表达式的值是后面项的值；逗号运算符优先级低于赋值运算符
    for (int i = int(n); i--;)
        if ((distance = spheres[i].intersect(r)) && distance < minDistance) {
            minDistance = distance;
            id = i;
        }
    // 所有内建赋值运算符都返回 *this
    return minDistance < inf;
}

Vec radiance (const Ray &r, int depth, unsigned short *Xi) { // 亮度函数，计算交点沿光线方向发射的亮度
    double distance;                               // distance to intersection
    int objID = 0;                               // objID of intersected object

    if (!intersect(r, distance, objID))
        return Vec(); // if miss, return black (0,0,0)

    const Sphere &obj = spheres[objID];        // the hit object

    if (depth > 10) return obj.emission;

    Vec position = r.origin + r.direction * distance;

    Vec normal = (position - obj.possion).normalize(); // 标准化
    Vec nl = normal.dot(r.direction) < 0 ? normal : normal * -1; // 判断方向，含方向的单位矢量
    Vec f = obj.color; // bsdf value

//    double p;
//    if (f.x > f.y && f.x > f.z)
//        p = f.x;
//    else if (f.y > f.z)
//        p = f.y;
//    else
//        p = f.z;

    // 得到f的最大的方向
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl

    if (++depth > 5)
            if (erand48(Xi) < p) f = f * (1 / p);
            else return obj.emission; //R.R.

    // 将面光源当点光源进行计算
    if (obj.refl == DIFF) {                  // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi);
        double r2 = erand48(Xi);
        double r2s = sqrt(r2);

        Vec w = nl;
        Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).normalize();
        Vec v = w % u;

        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize();
        return obj.emission + f.mult(radiance(Ray(position, d), depth, Xi));

    } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
        return obj.emission +
               f.mult(radiance(Ray(position, r.direction - normal * 2 * normal.dot(r.direction)), depth, Xi));

    Ray reflRay(position, r.direction - normal * 2 * normal.dot(r.direction));     // Ideal dielectric REFRACTION
    bool into = normal.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.direction.dot(nl), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)    // Total internal reflection
        return obj.emission + f.mult(radiance(reflRay, depth, Xi));
    Vec tdir = (r.direction * nnt - normal * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(normal));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);

    return obj.emission + f.mult(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette
                                              radiance(reflRay, depth, Xi) * RP :
                                              radiance(Ray(position, tdir), depth, Xi) *
                                              TP) :
                                 radiance(reflRay, depth, Xi) * Re + radiance(Ray(position, tdir), depth, Xi) * Tr);
}

int main (int argc, char *argv[]) {
    int width = 1024, height = 768;
    int samplesPerPixel = argc == 2 ? atoi(argv[1]) / 4 : 10; // # samples

    Ray camera(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize()); // camera pos, dir
    Vec cx = Vec(width * .5135 / height); // camera的左限
    Vec cy = (cx % camera.direction).normalize() * .5135; // camera的上限
    Vec r, *c = new Vec[width * height];

// 便捷的开启 OpenMP的多线程加速功能
#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y = 0; y < height; y++) {                       // Loop over image rows
        fprintf(stderr, "\rRendering (%direction spp) %5.2f%%", samplesPerPixel * 4, 100. * y / (height - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < width; x++)   // Loop cols

            // 2*2 MSAA抗锯齿
            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++)     // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()) {        // 2x2 subpixel cols

                    for (int s = 0; s < samplesPerPixel; s++) {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / height - .5) + camera.direction;

                        r = r + radiance(Ray(camera.origin + d * 140, d.normalize()), 0, Xi) * (1. / samplesPerPixel);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }

    // save by ppm
    FILE *f = fopen("image.ppm", "width");         // Write image to PPM file.
    fprintf(f, "P3\n%direction %direction\n%direction\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
        fprintf(f, "%direction %direction %direction ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}