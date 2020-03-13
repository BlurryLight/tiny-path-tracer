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
vec3 random_on_hemisphere();
vec3 random_on_sphere();

// thread-safe random generator
double drand_r(double min = 0.0, double max = 1.0);

vec3 reflect(const vec3 &v, const vec3 &n);
bool refract(const vec3 &v, const vec3 &n, float ni_over_nt, vec3 &refracted);
float schlick(float cosine, float ref_index);

inline float float_min(float a, float b) { return a < b ? a : b; }
inline float float_max(float a, float b) { return a > b ? a : b; }

vec3 color(const ray &r, hitable *world, hitable *light_shape, int depth,
           int max_depth);

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
class onb // ortho-normal bases
{
public:
  onb(const vec3 &normal) { build_from_w(normal); }
  inline vec3 operator[](int i) const { return axis_[i]; }
  vec3 u() const { return axis_[0]; }
  vec3 v() const { return axis_[1]; }
  vec3 w() const { return axis_[2]; }
  vec3 local(float x, float y, float z) const {
    return x * u() + y * v() + z * w();
  }
  vec3 local(const vec3 &a) { return a.x() * u() + a.y() * v() + a.z() * w(); }

  void build_from_w(const vec3 &normal);

  vec3 axis_[3];
};

class pdf // an abstract interface
{
public:
  virtual float value(const vec3 &direction) const = 0;
  virtual vec3 generate() const = 0;
};

class cosine_pdf : public pdf // lambertian
{
public:
  cosine_pdf(const vec3 &normal) { uvw_ = new onb(normal); }
  ~cosine_pdf() { delete uvw_; }
  virtual float value(const vec3 &direction) const override {
    float cosine = dot(unit_vector(direction), uvw_->w());
    if (cosine > 0) {
      return cosine;
    }
    return 0;
  }
  virtual vec3 generate() const override {
    return uvw_->local(random_on_hemisphere());
  }
  onb *uvw_;
};

class hitable_pdf : public pdf {
  // direction toward a hitable
public:
  hitable_pdf(hitable *ptr, const vec3 &origin) : ptr_(ptr), origin_(origin) {}
  virtual float value(const vec3 &direction) const override;
  virtual vec3 generate() const override;
  hitable *ptr_;
  vec3 origin_;
};

inline vec3 de_nan(const vec3 &vec) {
  vec3 tmp = vec;
  if (std::isnan(vec.x()))
    tmp[0] = 0;
  if (std::isnan(vec.y()))
    tmp[1] = 0;
  if (std::isnan(vec.z()))
    tmp[2] = 0;
  return tmp;
}

#endif
