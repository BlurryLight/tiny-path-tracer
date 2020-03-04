#include "hitable.h"
#include "aabb.h"
#include "utils.h"

bool lambertian::scatter(const ray &r_in, const hit_record &rec,
                         vec3 &attenuation, ray &scattered) const {
  vec3 target = rec.point + random_in_unit_sphere();
  scattered = ray(rec.point, target - rec.point, r_in.time());
  attenuation = albedo_;
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
static int box_x_compare(const void *a, const void *b) {
  AABB box_left, box_right;
  hitable *ah = *(hitable **)a;
  hitable *bh = *(hitable **)b;

  if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
    std::cerr << "no bounding box in bvh_node constructor\n";

  if (box_left.min().x() - box_right.min().x() < 0.0)
    return -1;
  else
    return 1;
}

static int box_y_compare(const void *a, const void *b) {
  AABB box_left, box_right;
  hitable *ah = *(hitable **)a;
  hitable *bh = *(hitable **)b;

  if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
    std::cerr << "no bounding box in bvh_node constructor\n";

  if (box_left.min().y() - box_right.min().y() < 0.0)
    return -1;
  else
    return 1;
}

static int box_z_compare(const void *a, const void *b) {
  AABB box_left, box_right;
  hitable *ah = *(hitable **)a;
  hitable *bh = *(hitable **)b;

  if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
    std::cerr << "no bounding box in bvh_node constructor\n";

  if (box_left.min().z() - box_right.min().z() < 0.0)
    return -1;
  else
    return 1;
}

bvh_node::bvh_node(hitable **l, int n, float time0, float time1) {
  int axis = int(drand_r(0, 3.0));
  if (axis == 0) {
    std::qsort(l, n, sizeof(hitable *), box_x_compare);
  } else if (axis == 1) {
    std::qsort(l, n, sizeof(hitable *), box_y_compare);
  } else {
    std::qsort(l, n, sizeof(hitable *), box_z_compare);
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
