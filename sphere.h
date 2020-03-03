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
  vec3 center_;
  float radius_;
  material *mat_ptr_;
};

class moving_sphere : public hitable {
public:
  moving_sphere() {}
  moving_sphere(vec3 center0, vec3 center1, float t0, float t1, float radius,
                material *m) {
    center0_ = center0;
    center1_ = center1;
    radius_ = radius;
    time0_ = t0;
    time1_ = t1;
    mat_ptr_ = m;
  }
  vec3 center(float time) const {
    return center0_ +
           (time - time0_) / (time1_ - time0_) * (center1_ - center0_);
  }
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;

  vec3 center0_, center1_;
  float radius_;
  material *mat_ptr_;
  float time0_, time1_;
};
bool sphere::hit(const ray &r, float t_min, float t_max,
                 hit_record &rec) const {
  vec3 oc = r.origin() - this->center_;
  float a = dot(r.direction(), r.direction());
  float b = 2.0 * dot(r.direction(), oc);
  float c = dot(oc, oc) - radius_ * radius_;

  float discriminant = b * b - 4 * a * c;

  if (discriminant > 0) {
    float temp = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.point = r.point_at_parameter(temp);
      rec.normal = (rec.point - center_) / radius_;
      rec.mat_ptr = mat_ptr_;
      return true;
    }
    temp = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.point = r.point_at_parameter(temp);
      rec.normal = (rec.point - center_) / radius_;
      rec.mat_ptr = mat_ptr_;
      return true;
    }
  }
  return false;
}

bool moving_sphere::hit(const ray &r, float t_min, float t_max,
                        hit_record &rec) const {
  vec3 oc = r.origin() - center(r.time());
  float a = dot(r.direction(), r.direction());
  float b = 2.0 * dot(r.direction(), oc);
  float c = dot(oc, oc) - radius_ * radius_;

  float discriminant = b * b - 4 * a * c;

  if (discriminant > 0) {
    float temp = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.point = r.point_at_parameter(temp);
      rec.normal = (rec.point - center(r.time())) / radius_;
      rec.mat_ptr = mat_ptr_;
      return true;
    }
    temp = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.point = r.point_at_parameter(temp);
      rec.normal = (rec.point - center(r.time())) / radius_;
      rec.mat_ptr = mat_ptr_;
      return true;
    }
  }
  return false;
}

#endif // SPHERE_H
