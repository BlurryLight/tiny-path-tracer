#ifndef CAMERA_H
#define CAMERA_H

#include "ray.h"
inline vec3 random_in_unit_disk() {
  vec3 p;
  do {
    p = 2.0 * vec3(drand_r(), drand_r(), 0) - vec3(1, 1, 0);
  } while (dot(p, p) >= 1.0);
  return p;
}
class camera {
public:
  // lookfrom is the origin
  // lookat is the point to look at
  // vup, the view up vector to project on the new plane when we incline it. We
  // can also tilt the plane
  camera(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov,
         float aspect) { // vfov is top to bottom in degrees, field of view on
                         // the vertical axis
    vec3 w, u, v;
    float theta = vfov * M_PI / 180; // convert to radiants
    float half_height = tan(theta / 2);
    float half_width = aspect * half_height;
    origin = lookfrom;
    w = unit_vector(lookfrom - lookat);
    u = unit_vector(cross(vup, w));
    v = cross(w, u);
    lower_left_corner = origin - u * half_width - v * half_height - w;
    horizontal = 2 * half_width * u;
    vertical = 2 * half_height * v;
  }
  ray get_ray(float s, float t) {
    return ray(origin,
               lower_left_corner + s * horizontal + t * vertical - origin);
  }

  vec3 origin;
  vec3 lower_left_corner;
  vec3 horizontal;
  vec3 vertical;
};

class camera_with_blur {
public:
  camera_with_blur(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov,
                   float aspect, float aperture, float focus_dist) {
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
