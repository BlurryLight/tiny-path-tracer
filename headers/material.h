#ifndef MATERIAL_H
#define MATERIAL_H

#include <ray.h>
#include <texture.h>
#include <vec3.h>

struct hit_record;

class material {
public:
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const = 0;
  virtual vec3 emitted(float u, float v, const vec3 &p) const {
    return vec3(0, 0, 0);
  }
};

class lambertian : public material {
public:
  lambertian(texture *albedo) : albedo_(albedo) {}
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override;

  //  vec3 albedo_; // reflection ratio
  texture *albedo_;
};

class metal : public material {
public:
  metal(const vec3 &albedo, float fuzz) : albedo_(albedo) {
    if (fuzz < 1 && fuzz >= 0)
      fuzz_ = fuzz;
    else
      fuzz_ = 1;
  }
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override;
  vec3 albedo_;
  float fuzz_;
  // fuzz_: 一个模糊系数，可以让金属的反射向量 = reflect + fuzz * random_unit
};
class dielectric : public material {
public:
  float ref_idx_;
  dielectric(float ri) : ref_idx_(ri) {}
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override;
};

class diffuse_light : public material {
public:
  diffuse_light(texture *a) : emit_(a) {}
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override {
    return false;
  }
  virtual vec3 emitted(float u, float v, const vec3 &p) const override {
    return emit_->value(u, v, p);
  }

  texture *emit_;
};

class isotropic : public material {
public:
  isotropic(texture *texture) : albedo_(texture) {}
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override;
  texture *albedo_;
};
#endif