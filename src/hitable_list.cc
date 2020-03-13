#include "hitable_list.h"
#include "aabb.h"
bool hitable_list::bounding_box(float t0, float t1, AABB &box) const {
  if (list_size_ < 1)
    return false;
  AABB temp_box;
  bool first_true = list_[0]->bounding_box(t0, t1, temp_box);
  if (!first_true) {
    // should never reach
    return false;
  } else {
    box = temp_box;
  }
  for (int i = 1; i < list_size_; i++) {
    if (list_[i]->bounding_box(t0, t1, temp_box)) {
      box = surrounding_box(box, temp_box);
    } else
      // should never reach
      return false;
  }
  return true;
}

float hitable_list::pdf_value(const vec3 &origin, const vec3 &direction) const {
  float weight = 1.0 / list_size_;
  float sum = 0;
  //权重和=1,概率密度积分也恒=1
  for (int i = 0; i < list_size_; i++)
    sum += weight * list_[i]->pdf_value(origin, direction);
  return sum;
}

vec3 hitable_list::random(const vec3 &origin) const {
  int index = int(drand_r() * list_size_);
  return list_[index]->random(origin);
}

bool hitable_list::hit(const ray &r, float t_min, float t_max,
                       hit_record &rec) const {
  hit_record temp_rec;
  bool hit_anything = false;
  double closest_so_far = t_max;
  for (int i = 0; i < list_size_; i++) {
    if (list_[i]->hit(r, t_min, closest_so_far, temp_rec)) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      rec = temp_rec;
    }
  }
  return hit_anything;
}
