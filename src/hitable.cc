#include "hitable.h"
#include "aabb.h"
#include "material.h"
#include "utils.h"
#include <array>
#include <hitable_list.h>
#include <limits>

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

constant_medium::constant_medium(hitable *boundary, float density,
                                 texture *texture)
    : boundary_(boundary), density_(density) {
  phase_funcion_ = new isotropic(texture);
}

bool constant_medium::hit(const ray &r, float t_min, float t_max,
                          hit_record &rec) const {

  hit_record rec1, rec2;
  if (boundary_->hit(r, std::numeric_limits<float>::lowest(),
                     std::numeric_limits<float>::max(),
                     rec1)) //可能从内向外命中
  {
    if (boundary_->hit(r, rec1.t + 0.0001, std::numeric_limits<float>::max(),
                       rec2)) {
      if (rec1.t < t_min)
        rec1.t = t_min;
      if (rec2.t > t_max)
        rec2.t = t_max;
      if (rec1.t >= rec2.t)
        return false;
      float distance_inside_boundary =
          (rec2.t - rec1.t) * r.direction().length();
      float hit_distance = (-1 / density_) * std::log(drand_r());
      if (hit_distance < distance_inside_boundary) {
        rec.t = rec1.t + hit_distance / r.direction().length();
        rec.point = r.point_at_parameter(rec.t);
        rec.normal = vec3(drand_r(), drand_r(), drand_r()); // arbirary vector
        rec.normal.normalize();
        rec.mat_ptr = phase_funcion_;
        return true;
      }
    }
  }
  return false;
}
