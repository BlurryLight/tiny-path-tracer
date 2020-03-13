#ifndef HITABLE_H
#define HITABLE_H

#include "aabb.h"
#include "ray.h"
#include "texture.h"
#include <algorithm>
// forward declaration
class AABB;
struct hit_record;
struct material;

struct hit_record {
  float t;
  float u = 0.0f;
  float v = 0.0f;
  vec3 point;
  vec3 normal;
  material *mat_ptr;
};

class hitable {
public:
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const = 0;
  virtual bool bounding_box(float t0, float t1, AABB &box) const = 0;
  virtual float pdf_value(const vec3 &origin, const vec3 &direction) const {
    return 0.0;
  }
  virtual vec3 random(const vec3 &origin) const { return vec3(1, 0, 0); }
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



class constant_medium : public hitable {
public:
  constant_medium(hitable *boundary, float density, texture *texture);
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  virtual bool bounding_box(float t0, float t1, AABB &box) const override {
    return boundary_->bounding_box(t0, t1, box);
  }
  hitable *boundary_;
  float density_;
  material *phase_funcion_;
};

#endif // HITABLE_H
