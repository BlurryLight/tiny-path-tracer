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
  virtual vec3 emitted(float u, float v, const vec3 &p) const {
    return vec3(0, 0, 0);
  }
};

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

class xy_rect : public hitable {
public:
  xy_rect() {}
  xy_rect(float x0, float x1, float y0, float y1, float k, material *mat)
      : x0_(x0), x1_(x1), y0_(y0), y1_(y1), k_(k), mat_ptr_(mat) {}
  virtual bool bounding_box(float t0, float t1, AABB &box) const override;
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  material *mat_ptr_;
  float x0_, x1_, y0_, y1_, k_; // z = k
};

class xz_rect : public hitable {
public:
  xz_rect() {}
  xz_rect(float x0, float x1, float z0, float z1, float k, material *mat)
      : x0_(x0), x1_(x1), z0_(z0), z1_(z1), k_(k), mat_ptr_(mat) {}
  virtual bool bounding_box(float t0, float t1, AABB &box) const override;
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  material *mat_ptr_;
  float x0_, x1_, z0_, z1_, k_; // y = k
};

class yz_rect : public hitable {
public:
  yz_rect() {}
  yz_rect(float y0, float y1, float z0, float z1, float k, material *mat)
      : y0_(y0), y1_(y1), z0_(z0), z1_(z1), k_(k), mat_ptr_(mat) {}
  virtual bool bounding_box(float t0, float t1, AABB &box) const override;
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  material *mat_ptr_;
  float y0_, y1_, z0_, z1_, k_; // y = k
};

class flip_normal : public hitable {
public:
  flip_normal(hitable *p) : ptr_(p) {}
  virtual bool bounding_box(float t0, float t1, AABB &box) const override {
    return ptr_->bounding_box(t0, t1, box);
  }
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override {
    if (ptr_->hit(r, t_min, t_max, rec)) {
      rec.normal = -rec.normal;
      return true;
    }
    return false;
  }

  hitable *ptr_;
};

class box : public hitable {
public:
  box() {}
  box(vec3 pmin, vec3 pmax, material *mat);
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  virtual bool bounding_box(float t0, float t1, AABB &box) const override {
    box = AABB(point_min_, point_max_);
    return true;
  }

  vec3 point_min_, point_max_;
  hitable *list_ptr_;
};

class translate : public hitable {
public:
  translate(hitable *p, const vec3 &offset) : ptr_(p), offset_(offset) {}
  virtual bool bounding_box(float t0, float t1, AABB &box) const override {
    if (ptr_->bounding_box(t0, t1, box)) {
      box = AABB(box.min() + offset_, box.max() + offset_);
      return true;
    }
    return false;
  }
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override {
    ray moved_r(r.origin() - offset_, r.direction(), r.time());
    if (ptr_->hit(moved_r, t_min, t_max, rec)) {
      rec.point += offset_;
      return true;
    }
    return false;
  }

  hitable *ptr_;
  vec3 offset_;
};

class rotate_y : public hitable {
public:
  rotate_y(hitable *p, float angle);
  virtual bool bounding_box(float t0, float t1, AABB &box) const override {
    box = box_;
    return has_box_;
  }
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;

  hitable *ptr_;
  float sin_theta_;
  float cos_theta_;
  bool has_box_;
  AABB box_;
};

#endif // HITABLE_H
