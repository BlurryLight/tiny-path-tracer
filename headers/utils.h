#ifndef UTILS_H
#define UTILS_H

#include "vec3.h"
struct hitable;
struct ray;
vec3 random_in_unit_disk();
vec3 random_in_unit_sphere();

// thread-safe random generator
double drand_r(double min = 0.0, double max = 1.0);

vec3 reflect(const vec3 &v, const vec3 &n);
bool refract(const vec3 &v, const vec3 &n, float ni_over_nt, vec3 &refracted);
float schlick(float cosine, float ref_index);

inline float float_min(float a, float b) { return a < b ? a : b; }
inline float float_max(float a, float b) { return a > b ? a : b; }

vec3 color(const ray &r, hitable *world, int depth, int max_depth);

hitable *random_scene();
#endif
