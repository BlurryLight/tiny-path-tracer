#ifndef UTILS_H
#define UTILS_H

#include "vec3.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <random>
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
hitable *light_spheres();
hitable *two_perlin_spheres();

hitable *sphere_cornell_box();
hitable *cornell_box();
hitable *cornell_box_smoke();
hitable *oneweek_final();

inline void get_uv_map(const vec3 &p, float &u, float &v) {
  u = std::atan2(p.z(), p.x()) / (2 * M_PI);
  v = std::asin(p.y()) / M_PI;
  u += 0.5f;
  v += 0.5f;
}

unsigned char *load_image_texture(std::string filename, int &width, int &height,
                                  int &channels);

#endif
