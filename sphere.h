#ifndef SPHERE_H
#define SPHERE_H

#include "hitable.h"

class sphere : public hitable {
public:
  sphere() {}
  sphere(vec3 center, float radius) : center_(center), radius_(radius) {}

  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;

  vec3 center_;
  float radius_;
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
      return true;
    }
    temp = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.point = r.point_at_parameter(temp);
      rec.normal = (rec.point - center_) / radius_;
      return true;
    }
  }
  return false;
}

#endif // SPHERE_H
