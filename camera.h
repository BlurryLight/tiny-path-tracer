#ifndef CAMERA_H
#define CAMERA_H

#include "ray.h"

class camera {
public:
  camera() {
    lower_left_corner_ = vec3(-2.0, -1.0, -1.0);
    horizontal_ = vec3(4.0, 0.0, 0.0);
    vertical_ = vec3(0.0, 2.0, 0.0);
    origin_ = vec3(0.0, 0.0, 0.0);
  }

  ray get_ray(float u, float v) {
    return ray(origin_,
               lower_left_corner_ + u * horizontal_ + v * vertical_ - origin_);
  }

  vec3 origin_;
  vec3 lower_left_corner_;
  vec3 vertical_;
  vec3 horizontal_;
};

#endif // CAMERA_H
