#include "rect_box.h"
#include "hitable_list.h"
bool xy_rect::bounding_box(float t0, float t1, AABB &box) const {
  box = AABB(vec3(x0_, y0_, k_ - 0.0001), vec3(x1_, y1_, k_ + 0.0001));
  return true;
}

bool xy_rect::hit(const ray &r, float t_min, float t_max,
                  hit_record &rec) const {
  float t = (k_ - r.origin().z()) / r.direction().z();
  if (t > t_max || t < t_min)
    return false;
  float x = r.origin().x() + t * r.direction().x();
  float y = r.origin().y() + t * r.direction().y();
  if (x < x0_ || x > x1_ || y < y0_ || y > y1_)
    return false;
  rec.u = (x - x0_) / (x1_ - x0_);
  rec.v = (y - y0_) / (y1_ - y0_);
  rec.t = t;
  rec.mat_ptr = mat_ptr_;
  rec.point = r.point_at_parameter(t);
  rec.normal = vec3(0, 0, 1);
  return true;
}

float xz_rect::pdf_value(const vec3 &origin, const vec3 &direction) const {
  hit_record rec;
  if (this->hit(ray(origin, direction), 0.0001,
                std::numeric_limits<float>::max(), rec)) {
    float rec_area = std::abs((x1_ - x0_) * (z1_ - z0_));
    float distance_squared =
        (rec.t * direction).squared_length(); // is that right?
    float cosine = std::abs(dot(direction, rec.normal) / direction.length());
    return distance_squared / (cosine * rec_area);
  }
  return 0.0f;
}

vec3 xz_rect::random(const vec3 &origin) const {
  vec3 random_point =
      vec3(x0_ + drand_r() * (x1_ - x0_), k_, z0_ + drand_r() * (z1_ - z0_));
  return random_point - origin;
}

bool xz_rect::bounding_box(float t0, float t1, AABB &box) const {

  box = AABB(vec3(x0_, k_ - 0.0001, z0_), vec3(x1_, k_ + 0.0001, z1_));
  return true;
}

bool xz_rect::hit(const ray &r, float t_min, float t_max,
                  hit_record &rec) const {

  float t = (k_ - r.origin().y()) / r.direction().y();
  if (t > t_max || t < t_min)
    return false;
  float x = r.origin().x() + t * r.direction().x();
  float z = r.origin().z() + t * r.direction().z();
  if (x < x0_ || x > x1_ || z < z0_ || z > z1_)
    return false;
  rec.u = (x - x0_) / (x1_ - x0_);
  rec.v = (z - z0_) / (z1_ - z0_);
  rec.t = t;
  rec.mat_ptr = mat_ptr_;
  rec.point = r.point_at_parameter(t);
  rec.normal = vec3(0, 1, 0);
  return true;
}

bool yz_rect::bounding_box(float t0, float t1, AABB &box) const {
  box = AABB(vec3(k_ - 0.0001, y0_, z0_), vec3(k_ + 0.0001, y1_, z1_));
  return true;
}

bool yz_rect::hit(const ray &r, float t_min, float t_max,
                  hit_record &rec) const {
  float t = (k_ - r.origin().x()) / r.direction().x();
  if (t > t_max || t < t_min)
    return false;
  float y = r.origin().y() + t * r.direction().y();
  float z = r.origin().z() + t * r.direction().z();
  if (y < y0_ || y > y1_ || z < z0_ || z > z1_)
    return false;
  rec.u = (y - y0_) / (y1_ - y0_);
  rec.v = (z - z0_) / (z1_ - z0_);
  rec.t = t;
  rec.mat_ptr = mat_ptr_;
  rec.point = r.point_at_parameter(t);
  rec.normal = vec3(1, 0, 0);
  return true;
}

box::box(vec3 pmin, vec3 pmax, material *mat) {
  point_min_ = pmin;
  point_max_ = pmax;
  hitable **list = new hitable *[6];
  // normals outword from the box

  // xy group
  list[0] = new xy_rect(pmin.x(), pmax.x(), pmin.y(), pmax.y(), pmax.z(), mat);
  list[1] = new flip_normal(
      new xy_rect(pmin.x(), pmax.x(), pmin.y(), pmax.y(), pmin.z(), mat));

  // xz group
  list[2] = new xz_rect(pmin.x(), pmax.x(), pmin.z(), pmax.z(), pmax.y(), mat);
  list[3] = new flip_normal(
      new xz_rect(pmin.x(), pmax.x(), pmin.z(), pmax.z(), pmin.y(), mat));

  // yz group
  list[4] = new yz_rect(pmin.y(), pmax.y(), pmin.z(), pmax.z(), pmax.x(), mat);
  list[5] = new flip_normal(
      new yz_rect(pmin.y(), pmax.y(), pmin.z(), pmax.z(), pmin.x(), mat));

  this->list_ptr_ = new hitable_list(list, 6);
}

bool box::hit(const ray &r, float t_min, float t_max, hit_record &rec) const {
  return list_ptr_->hit(r, t_min, t_max, rec);
}

rotate_y::rotate_y(hitable *p, float angle) : ptr_(p) {
  float radians = angle / 180.0f * M_PI;
  sin_theta_ = std::sin(radians);
  cos_theta_ = std::cos(radians);
  has_box_ = ptr_->bounding_box(0, 0, box_); // don't support moving-shpere
  auto fmax = std::numeric_limits<float>::max();
  vec3 min(fmax, fmax, fmax);
  vec3 max(-fmax, -fmax, -fmax);
  // find new bbox(new vec3 min . vec3 max)
  //   naive method : original book version
  //   for (int i = 0; i < 2; i++) {
  //     for (int j = 0; j < 2; j++) {

  //       for (int k = 0; k < 2; k++) {
  //         float x = i * box_.max().x() + (1 - i) * box_.min().x();
  //         float y = j * box_.max().y() + (1 - j) * box_.min().y();
  //         float z = k * box_.max().z() + (1 - k) * box_.min().z();
  //         float new_x = cos_theta_ * x + sin_theta_ * z;
  //         float new_z = -sin_theta_ * x + cos_theta_ * z;
  //         vec3 tester(new_x, y, new_z);
  //         for (int m = 0; m < 3; m++) {
  //           if (tester[m] > max[m])
  //             max[m] = tester[m];
  //           if (tester[m] < min[m])
  //             min[m] = tester[m];
  //         }
  //       }
  //     }
  //   }
  // effecitve version
  // see http://dev.theomader.com/transform-bounding-boxes/
  vec3 first_column{cos_theta_, 0, -sin_theta_};
  vec3 second_column{0, 1, 0};
  vec3 third_column{sin_theta_, 0, cos_theta_};

  vec3 xa = first_column * box_.max().x();
  vec3 xb = first_column * box_.min().x();

  vec3 ya = second_column * box_.max().y();
  vec3 yb = second_column * box_.min().y();

  vec3 za = third_column * box_.max().z();
  vec3 zb = third_column * box_.min().z();

  min = vec_min(xa, xb) + vec_min(ya, yb) + vec_min(za, zb);
  max = vec_max(xa, xb) + vec_max(ya, yb) + vec_max(za, zb);

  box_ = AABB(min, max);
}

bool rotate_y::hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const {
  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin[0] = cos_theta_ * r.origin()[0] - sin_theta_ * r.origin()[2];
  origin[2] = sin_theta_ * r.origin()[0] + cos_theta_ * r.origin()[2];

  direction[0] = cos_theta_ * r.direction()[0] - sin_theta_ * r.direction()[2];
  direction[2] = sin_theta_ * r.direction()[0] + cos_theta_ * r.direction()[2];
  ray rotated_r(origin, direction, r.time());
  if (ptr_->hit(rotated_r, t_min, t_max, rec)) {
    vec3 p = rec.point;
    vec3 normal = rec.normal; // if scale, the normal will has a more
                              // complicated computation
    p[0] = cos_theta_ * rec.point[0] + sin_theta_ * rec.point[2];
    p[2] = -sin_theta_ * rec.point[0] + cos_theta_ * rec.point[2];

    normal[0] = cos_theta_ * rec.normal[0] + sin_theta_ * rec.normal[2];
    normal[2] = -sin_theta_ * rec.normal[0] + cos_theta_ * rec.normal[2];
    rec.point = p;
    rec.normal = normal;
    return true;
  }
  return false;
}
