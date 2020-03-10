#ifndef PERLIN_NOISE_H
#define PERLIN_NOISE_H
#include <array>
#include <random>
#include <vec3.h>

class perlin_noise {
public:
  perlin_noise();
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
