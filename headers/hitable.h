#ifndef HITABLE_H
#define HITABLE_H

#include "aabb.h"
#include "ray.h"
#include "texture.h"
#include <algorithm>
// forward declaration
class AABB;
struct hit_record;

class material {
public:
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const = 0;
};

struct hit_record {
  float t;
  vec3 point;
  vec3 normal;
  material *mat_ptr;
};

class hitable {
public:
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const = 0;
  virtual bool bounding_box(float t0, float t1, AABB &box) const = 0;
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


class bvh_node : public hitable // binary tree
{
public:
  bvh_node() {}
  bvh_node(hitable **l, int n, float time0, float time1);
  // hitable** is a pointer to an array of hitbale
  virtual bool bounding_box(float t0, float t1, AABB &box) const override;
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  hitable *left_;  // hitable object( BVH_nodes,spheres and others)
  hitable *right_; //
  AABB box_;
};


#endif // HITABLE_H
