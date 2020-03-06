#include "sphere.h"

moving_sphere::moving_sphere(vec3 center0, vec3 center1, float t0, float t1,
                             float radius, material *m) {
  center0_ = center0;
  center1_ = center1;
  radius_ = radius;
  time0_ = t0;
  time1_ = t1;
  mat_ptr_ = m;
}

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
      get_uv_map((rec.point - center_) / radius_, rec.u, rec.v);
      rec.normal = (rec.point - center_) / radius_;
      rec.mat_ptr = mat_ptr_;
      return true;
    }
    temp = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.point = r.point_at_parameter(temp);
      get_uv_map((rec.point - center_) / radius_, rec.u, rec.v);
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
      get_uv_map((rec.point - center0_) / radius_, rec.u, rec.v);
      rec.mat_ptr = mat_ptr_;
      return true;
    }
    temp = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.point = r.point_at_parameter(temp);
      rec.normal = (rec.point - center(r.time())) / radius_;
      get_uv_map((rec.point - center0_) / radius_, rec.u, rec.v);
      rec.mat_ptr = mat_ptr_;
      return true;
    }
  }
  return false;
}

bool moving_sphere::bounding_box(float t0, float t1, AABB &box) const {

  AABB box0{center(t0) - vec3(radius_, radius_, radius_),
            center(t0) + vec3(radius_, radius_, radius_)};
  AABB box1{center(t1) - vec3(radius_, radius_, radius_),
            center(t1) + vec3(radius_, radius_, radius_)};
  box = surrounding_box(box0, box1);
  return true;
}

bool sphere::bounding_box(float t0, float t1, AABB &box) const {
  box = AABB(center_ - vec3(radius_, radius_, radius_),
             center_ + vec3(radius_, radius_, radius_));
  return true;
}
