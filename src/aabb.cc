#include "aabb.h"

bool AABB::hit(const ray &r, float tmin, float tmax) const {
  for (int i = 0; i < 3; i++) {
    float inv_D = 1.0f / r.direction()[i];
    float t0 = (min_[i] - r.origin()[i]) * inv_D;
    float t1 = (max_[i] - r.origin()[i]) * inv_D;
    if (inv_D < 0.0f) {
      std::swap(t0, t1);
    }
    tmin = t0 > tmin ? t0 : tmin;
    tmax = t1 < tmax ? t1 : tmax;
    if (tmin > tmax) {
      // no overlap
      return false;
    }
  }
  return true;
}

AABB surrounding_box(AABB box0, AABB box1) {
  vec3 small{float_min(box0.min().x(), box1.min().x()),
             float_min(box0.min().y(), box1.min().y()),
             float_min(box0.min().z(), box1.min().z())};

  vec3 big{float_max(box0.max().x(), box1.max().x()),
           float_max(box0.max().y(), box1.max().y()),
           float_max(box0.max().z(), box1.max().z())};
  return AABB(small, big);
}
