#ifndef RAY_H
#define RAY_H
#include "vec3.h"
#include <random>
class ray {
public:
  ray();
  ray(const vec3 &pos, const vec3 &direction, float time = 0.0f)
      : pos_(pos), direction_(direction), time_(time) {}
  vec3 origin() const { return pos_; }
  vec3 direction() const { return direction_; }
  float time() const { return time_; }
  vec3 point_at_parameter(float t) const { return pos_ + t * direction_; }

  vec3 pos_;
  vec3 direction_;
  float time_;
};

#endif // RAY_H
