#include "hitable.h"
#include "aabb.h"
#include "utils.h"
#include <array>
#include <hitable_list.h>

bool lambertian::scatter(const ray &r_in, const hit_record &rec,
                         vec3 &attenuation, ray &scattered) const {
  vec3 target = rec.point + random_in_unit_sphere();
  scattered = ray(rec.point, target - rec.point, r_in.time());
  attenuation = albedo_->value(rec.u, rec.v, rec.point);
  return true;
}

bool metal::scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
                    ray &scattered) const {
  vec3 reflected = reflect(unit_vector(r_in.direction_), rec.normal);
  scattered =
      ray(rec.point, reflected + fuzz_ * random_in_unit_sphere(), r_in.time());
  attenuation = albedo_;
  return (dot(scattered.direction(), rec.normal) > 0);
}

bool dielectric::scatter(const ray &r_in, const hit_record &rec,
                         vec3 &attenuation, ray &scattered) const {
  vec3 outward_normal;
  vec3 reflected = reflect(r_in.direction(), rec.normal);
  float ni_over_nt;
  attenuation = vec3(1.0, 1.0, 1.0);

  vec3 refracted;
  float reflect_prob;
  float cosine;

  if (dot(r_in.direction(), rec.normal) > 0) {
    //从内部射往外面
    outward_normal = -rec.normal;
    ni_over_nt = ref_idx_;
    // why there is a ref_idx?
    // https://zhuanlan.zhihu.com/p/47991519
    // 原因是schlick近似只能取表面外的角
    // 由内部向外，n1 = ref_ind, n2 = 1.0(空气)，cosine(theta_1)已知
    // 已知 n1 sin theta_1 = n2 sin theta_2
    // 做个近似，令sin theta_1 = cos_theta_1
    // 就可以推出来下面的式子了
    // 我更觉得这像一个错误
    // wrong version:
    // cosine = ref_idx_ * dot(r_in.direction(), rec.normal) /
    //          r_in.direction().length();
    // corret version:
    cosine = dot(r_in.direction(), rec.normal) / r_in.direction().length();
    cosine = std::sqrt(1 - ref_idx_ * ref_idx_ * (1 - cosine * cosine));

  } else {

    outward_normal = rec.normal;
    ni_over_nt = 1.0 / ref_idx_; // 1.0空气的折射率
    cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
  }
  if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
    reflect_prob = schlick(cosine, ref_idx_);
  } else {
    reflect_prob = 1.0;
    //全反射
  }
  //反射概率
  //因为我们要采样多次，所以这里可以用随机数来模拟概率(蒙特卡罗)
  if (drand_r() < reflect_prob) {
    scattered = ray(rec.point, reflected, r_in.time());
  } else {
    scattered = ray(rec.point, refracted, r_in.time());
  }
  return true;
}

// BVH compare operator
struct box_axis_compare {
  // 0 x_axis_compare
  // 1 y
  // 2 z
  //<0 || > 2 ub
  box_axis_compare(int axis) : axis_(axis) {}
  bool operator()(const hitable *ah, const hitable *bh) {
    AABB box_left, box_right;
    auto a = ah->bounding_box(0, 0, box_left);
    auto b = bh->bounding_box(0, 0, box_right);
    if (!a || !b) {
      std::cout << "No bounding box in BVH node constructor" << std::endl;
      exit(-1);
    }
    if (std::isgreaterequal(box_left.min()[axis_], box_right.min()[axis_])) {
      return false;
    }
    return true;
  }
  int axis_ = 0;
};

bvh_node::bvh_node(hitable **l, int n, float time0, float time1) {
  int axis = int(drand_r(0, 3.0));
  if (axis == 0) {
    std::sort(l, l + n, box_axis_compare(0));
  } else if (axis == 1) {
    std::sort(l, l + n, box_axis_compare(1));
  } else {
    std::sort(l, l + n, box_axis_compare(2));
  }
  if (n == 1) {
    left_ = right_ = l[0];
  } else if (n == 2) {
    left_ = l[0];
    right_ = l[1];
  } else {
    left_ = new bvh_node(l, n / 2, time0, time1);
    right_ = new bvh_node(l + n / 2, n - n / 2, time0, time1);
  }
  AABB box_left, box_right;
  if (!left_->bounding_box(time0, time1, box_left) ||
      !right_->bounding_box(time0, time1, box_right)) {
    std::cout << "no bounding box in bvh constructor" << std::endl;
  }
  box_ = surrounding_box(box_left, box_right);
}

bool bvh_node::bounding_box(float t0, float t1, AABB &box) const {
  box = box_;
  return true;
}

bool bvh_node::hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const {
  if (box_.hit(r, t_min, t_max)) {
    hit_record left_rec, right_rec;
    bool hit_left = left_->hit(r, t_min, t_max, left_rec);
    bool hit_right = right_->hit(r, t_min, t_max, right_rec);
    if (hit_left && hit_right) // hit both
    {
      if (left_rec.t < right_rec.t) // left objects are more closer to camera
      {
        rec = left_rec;
      } else {
        rec = right_rec;
      }
      return true;
    } else if (hit_left) {
      rec = left_rec;
      return true;
    } else if (hit_right) {
      rec = right_rec;
      return true;
    } else {
      return false;
    }
  }
  // doesn't hit box
  return false;
}

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
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {

      for (int k = 0; k < 2; k++) {
        float x = i * box_.max().x() + (1 - i) * box_.min().x();
        float y = j * box_.max().y() + (1 - j) * box_.min().y();
        float z = k * box_.max().z() + (1 - k) * box_.min().z();
        float new_x = cos_theta_ * x + sin_theta_ * z;
        float new_z = -sin_theta_ * x + cos_theta_ * z;
        vec3 tester(new_x, y, new_z);
        for (int m = 0; m < 3; m++) {
          if (tester[m] > max[m])
            max[m] = tester[m];
          if (tester[m] < min[m])
            min[m] = tester[m];
        }
      }
    }
  }
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
