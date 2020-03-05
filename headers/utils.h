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
    std::generate(random_vec3_.begin(), random_vec3_.end(), []() -> vec3 {
      float f1 = static_cast<float>(2 * drand_r() - 1);
      float f2 = static_cast<float>(2 * drand_r() - 1);
      float f3 = static_cast<float>(2 * drand_r() - 1);
      return unit_vector(vec3(f1, f2, f3));
    });
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
  float turb(const vec3 &p, int depth = 5) const;

  static std::array<vec3, 256> random_vec3_;
  static std::array<int, 256> permute_x_; // std::shuffle
  static std::array<int, 256> permute_y_; // std::shuffle
  static std::array<int, 256> permute_z_; // std::shuffle
private:
  float perlin_interpolate(vec3 vertex[2][2][2], float u, float v,
                           float w) const;
};

#endif
