#ifndef CAMERA_H
#define CAMERA_H

#include "ray.h"
#include "utils.h"
class camera_with_blur {
public:
  camera_with_blur(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov,
                   float aspect, float aperture, float focus_dist, float t0,
                   float t1);

  ray get_ray(float s, float t);

  vec3 origin_;
  vec3 lower_left_corner_;
  vec3 vertical_;
  vec3 horizontal_;
  vec3 u_, v_, w_;
  float lens_radius_;
  float time0, time1; // aperture open time
};

using camera = camera_with_blur;

#endif // CAMERA_H
