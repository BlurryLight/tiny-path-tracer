#include "material.h"
#include "hitable.h"
bool lambertian::scatter(const ray &r_in, const hit_record &hit_rec,
                         scatter_record &scatter_rec) const {
  scatter_rec.is_specular = false;
  scatter_rec.attenuation = albedo_->value(hit_rec.u, hit_rec.v, hit_rec.point);
  scatter_rec.pdf_ptr.reset(new cosine_pdf(hit_rec.normal));
  return true;
}

float lambertian::scattering_pdf(const ray &r_in, const hit_record &rec,
                                 const ray &scattered) const {
  float cosine = dot(rec.normal, unit_vector(scattered.direction()));
  if (cosine < 0)
    return 0;
  return cosine / M_PI;
}

bool dielectric::scatter(const ray &r_in, const hit_record &rec,
                         scatter_record &scatter_rec) const {
  vec3 outward_normal;
  vec3 reflected = reflect(r_in.direction(), rec.normal);
  float ni_over_nt;
  scatter_rec.is_specular = true;
  scatter_rec.attenuation = vec3(1.0, 1.0, 1.0);

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
    scatter_rec.specular_ray = ray(rec.point, reflected, r_in.time());
  } else {
    scatter_rec.specular_ray = ray(rec.point, refracted, r_in.time());
  }
  return true;
}

bool isotropic::scatter(const ray &r_in, const hit_record &rec,
                        vec3 &attenuation, ray &scattered) const {
  scattered = ray(rec.point, random_in_unit_sphere());
  attenuation = albedo_->value(rec.u, rec.v, rec.point);
  return true;
}

vec3 diffuse_light::emitted(const ray &r_in, const hit_record &rec, float u,
                            float v, const vec3 &p) const {
  if (dot(rec.normal, r_in.direction()) < 0)
    return emit_->value(u, v, p);
  else {
    return vec3(0, 0, 0);
  }
}

bool metal::scatter(const ray &r_in, const hit_record &rec,
                    scatter_record &scatter_rec) const {

  vec3 reflected = reflect(unit_vector(r_in.direction_), rec.normal);
  scatter_rec.specular_ray =
      ray(rec.point, reflected + fuzz_ * random_in_unit_sphere(), r_in.time());
  scatter_rec.attenuation = albedo_;
  scatter_rec.is_specular = true;
  scatter_rec.pdf_ptr = nullptr;
  return (dot(scatter_rec.specular_ray.direction(), rec.normal) > 0);
}
