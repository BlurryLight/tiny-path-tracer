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
hitable *two_spheres();
hitable *two_perlin_spheres();

class perlin_noise {
public:
  perlin_noise() {
    std::generate(random_float_.begin(), random_float_.end(),
                  []() { return drand_r(); });
    for (int i = 0; i != permute_x_.size(); i++) {
      permute_x_[i] = i;
      permute_y_[i] = i;
      permute_z_[i] = i;
    }
    auto seed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    auto engine = std::default_random_engine(seed);
    std::shuffle(permute_x_.begin(), permute_x_.end(), engine);
    std::shuffle(permute_y_.begin(), permute_y_.end(), engine);
    std::shuffle(permute_z_.begin(), permute_z_.end(), engine);
  }
  float noise(const vec3 &p) const;

  static std::array<float, 256> random_float_;
  static std::array<int, 256> permute_x_; // std::shuffle
  static std::array<int, 256> permute_y_; // std::shuffle
  static std::array<int, 256> permute_z_; // std::shuffle
private:
  float trilinear_interpolate(float vertex[2][2][2], float u, float v,
                              float w) const;
  float perlin_interpolate(float vertex[2][2][2], float u, float v,
                           float w) const;
};

#endif
