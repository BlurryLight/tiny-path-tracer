#include "camera.h"
camera_with_blur::camera_with_blur(vec3 lookfrom, vec3 lookat, vec3 vup,
                                   float vfov, float aspect, float aperture,
                                   float focus_dist, float t0, float t1) {
  // aperture = 0.0, edge-blur disappears
  // time0 = time1 = 0.0, motion-blur off
  time0 = t0;
  time1 = t1;
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

ray camera_with_blur::get_ray(float s, float t) {
  vec3 rd = lens_radius_ * random_in_unit_disk();
  vec3 offset = u_ * rd.x() + v_ * rd.y();
  float time = time0 + drand_r() * (time1 - time0);
  return ray(origin_ + offset,
             lower_left_corner_ + s * horizontal_ + t * vertical_ - origin_ -
                 offset,
             time);
}
