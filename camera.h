#ifndef CAMERA_H
#define CAMERA_H

#include "ray.h"
inline vec3 random_in_unit_disk() {
  vec3 p;
  do {
    p = 2.0 * vec3(drand48(), drand48(), 0) - vec3(1, 1, 0);
  } while (dot(p, p) >= 1.0);
  return p;
}
class camera {
public:
  camera(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov, float aspect,
         float aperture, float focus_dist) {
    lens_radius_ = aperture / 2;
    float theta = vfov * M_PI / 180;
    float half_height = tan(theta / 2);
    float half_width = aspect * half_height;
    origin_ = lookfrom;
    w_ = unit_vector(lookfrom - lookat);
    u_ = unit_vector(cross(vup, w_));
    v_ = unit_vector(cross(w_, u_));
    lower_left_corner_ = origin_ - half_width * focus_dist * u_ -
                         half_height * focus_dist * v_ - w_ * focus_dist;
    horizontal_ = 2 * half_width * u_ * focus_dist;
    vertical_ = 2 * half_height * v_ * focus_dist;
  }

  ray get_ray(float s, float t) {
    vec3 rd = lens_radius_ * random_in_unit_disk();
    vec3 offset = u_ * rd.x() + v_ * rd.y();
    return ray(origin_ + offset, lower_left_corner_ + s * horizontal_ +
                                     t * vertical_ - origin_ - offset);
  }

  vec3 origin_;
  vec3 lower_left_corner_;
  vec3 vertical_;
  vec3 horizontal_;
  vec3 u_, v_, w_;
  float lens_radius_;
};

#endif // CAMERA_H
