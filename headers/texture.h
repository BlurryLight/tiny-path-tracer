#ifndef TEXTURE_H
#define TEXTURE_H
#include "utils.h"
#include "vec3.h"
#include <iostream>

struct perlin_noise;

class texture {
public:
  virtual vec3 value(float u, float v, const vec3 &p) const = 0;
};
// pure color

class constant_texture : public texture {
public:
  constant_texture() {}
  constant_texture(vec3 c) : color_(c) {}

  virtual vec3 value(float u, float v, const vec3 &p) const override {
    return color_;
  }
  vec3 color_;
};
class checker_texture : public texture {
public:
  checker_texture() {}
  checker_texture(texture *t0, texture *t1) : odd_(t0), even_(t1) {}

  virtual vec3 value(float u, float v, const vec3 &p) const override;
  texture *odd_;
  texture *even_;
};

class perlin_noise_texture : public texture {
public:
  perlin_noise_texture() {}
  virtual vec3 value(float u, float v, const vec3 &p) const override;
  perlin_noise noise_;
};

#endif
