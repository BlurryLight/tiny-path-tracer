#ifndef AABB_H
#define AABB_H

#include "ray.h"
#include "utils.h"
#include "vec3.h"
class AABB {
public:
  AABB() {}
  AABB(const vec3 &a, const vec3 &b) {
    min_ = a;
    max_ = b;
  }
  vec3 min() const { return min_; }
  vec3 max() const { return max_; }
  bool hit(const ray &r, float tmin, float tmax) const;
  vec3 min_;
  vec3 max_;
};

AABB surrounding_box(AABB box0, AABB box1);
#endif // AABB_H
