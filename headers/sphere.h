#ifndef SPHERE_H
#define SPHERE_H

#include "hitable.h"

class sphere : public hitable {
public:
  sphere() {}
  sphere(vec3 center, float radius, material *m)
      : center_(center), radius_(radius), mat_ptr_(m) {}

  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  virtual bool bounding_box(float t0, float t1, AABB &box) const override;
  virtual float pdf_value(const vec3 &origin,
                          const vec3 &direction) const override;
  virtual vec3 random(const vec3 &origin) const override;
  vec3 center_;
  float radius_;
  material *mat_ptr_;
};

class moving_sphere : public hitable {
public:
  moving_sphere() {}
  moving_sphere(vec3 center0, vec3 center1, float t0, float t1, float radius,
                material *m);
  vec3 center(float time) const {
    return center0_ +
           (time - time0_) / (time1_ - time0_) * (center1_ - center0_);
  }
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  virtual bool bounding_box(float t0, float t1, AABB &box) const override;

  vec3 center0_, center1_;
  float radius_;
  material *mat_ptr_;
  float time0_, time1_;
};

#endif // SPHERE_H
