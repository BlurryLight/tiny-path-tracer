#ifndef HITABLE_H
#define HITABLE_H

#include "ray.h"

inline vec3 random_in_unit_sphere() {
  vec3 p;
  do {
    p = 2.0 * vec3(drand48(), drand48(), drand48()) - vec3(1.0, 1.0, 1.0);
  } while (p.squared_length() >= 1.0);
  return p;
}
struct hit_record;
class material {
public:
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const = 0;
};

struct hit_record {
  float t;
  vec3 point;
  vec3 normal;
  material *mat_ptr;
};
class hitable {
public:
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const = 0;
};

class lambertian : public material {
public:
  lambertian(const vec3 &albedo) : albedo_(albedo) {}
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override {
    vec3 target = rec.point + random_in_unit_sphere();
    scattered = ray(rec.point, target - rec.point);
    attenuation = albedo_;
    return true;
  }

  vec3 albedo_; // reflection ratio
};

inline vec3 reflect(const vec3 &v, const vec3 &n) // vecin and normal
{
  return v - 2 * dot(v, n) * n;
}
class metal : public material {
public:
  metal(const vec3 &albedo, float fuzz) : albedo_(albedo) {
    if (fuzz < 1 && fuzz >= 0)
      fuzz_ = fuzz;
    else
      fuzz_ = 1;
    //一个模糊系数，可以让金属的反射向量 = reflect + fuzz * random_unit
  }
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override {
    vec3 reflected = reflect(unit_vector(r_in.direction_), rec.normal);
    scattered = ray(rec.point, reflected + fuzz_ * random_in_unit_sphere());
    attenuation = albedo_;
    return (dot(scattered.direction(), rec.normal) > 0);
  }
  vec3 albedo_;
  float fuzz_;
};
inline bool refract(const vec3 &v, const vec3 &n, float ni_over_nt,
                    vec3 &refracted) {
  vec3 unit_v = unit_vector(v);
  float dt = dot(unit_v, n);
  float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if (discriminant > 0) {
    refracted = ni_over_nt * (unit_v - n * dt) - n * sqrt(discriminant);
    return true;
  } else {
    return false;
  }
}
inline float schlick(float cosine, float ref_index) {
  // cosine 是ray_in和normal的dot结果,ref_index 折射系数,得到的结果是反射率
  float r0 = (1 - ref_index) / (1 + ref_index);
  r0 = r0 * r0;
  return r0 + (1 - r0) * pow((1 - cosine), 5);
}
class dielectric : public material {
public:
  float ref_idx_;
  dielectric(float ri) : ref_idx_(ri) {}
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       vec3 &attenuation, ray &scattered) const override {
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
    if (drand48() < reflect_prob) {
      scattered = ray(rec.point, reflected);
    } else {
      scattered = ray(rec.point, refracted);
    }
    return true;
  }
};

#endif // HITABLE_H
