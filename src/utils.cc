#include "utils.h"
#include "hitable.h"
#include "hitable_list.h"
#include "ray.h"
#include <random>
#include <sphere.h>
vec3 random_in_unit_disk() {
  vec3 p;
  do {
    p = 2.0 * vec3(drand_r(), drand_r(), 0) - vec3(1, 1, 0);
  } while (dot(p, p) >= 1.0);
  return p;
}

vec3 random_in_unit_sphere() {
  vec3 p;
  do {
    p = 2.0 * vec3(drand_r(), drand_r(), drand_r()) - vec3(1.0, 1.0, 1.0);
  } while (p.squared_length() >= 1.0);
  return p;
}
double drand_r(double min, double max) {
  static thread_local std::mt19937 generator;
  std::uniform_real_distribution<double> dis(min, max); //[min,max)
  return dis(generator);
}

bool refract(const vec3 &v, const vec3 &n, float ni_over_nt, vec3 &refracted) {
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

float schlick(float cosine, float ref_index) {
  // cosine 是ray_in和normal的dot结果,ref_index 折射系数,得到的结果是反射率
  float r0 = (1 - ref_index) / (1 + ref_index);
  r0 = r0 * r0;
  return r0 + (1 - r0) * pow((1 - cosine), 5);
}

vec3 reflect(const vec3 &v, const vec3 &n) // vecin and normal
{
  return v - 2 * dot(v, n) * n;
}

vec3 color(const ray &r, hitable *world, int depth, int max_depth) {
  hit_record rec;
  if (world->hit(r, 0.001, MAXFLOAT, rec)) {
    ray scattered;
    vec3 attenuation;
    if (depth < max_depth &&
        rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
      return attenuation * color(scattered, world, depth + 1, max_depth);
    } else {
      return vec3(0, 0, 0);
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = (unit_direction.y() + 1.0) * 0.5; // clamp (-1,1) to (0,1)
    return (1 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
  }
  // linear interperation
  // blend white with blue
}

hitable *random_scene() {
  int n = 500;
  hitable **list = new hitable *[n + 1];
  texture *checker =
      new checker_texture(new constant_texture({0.1, 0.1, 0.1}), // white
                          new constant_texture({0.9, 0.9, 0.9})  // black
      );
  list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(checker));
  int i = 1;
  for (int a = -11; a < 11; a++) {
    for (int b = -11; b < 11; b++) {
      float choose_mat = drand_r();
      vec3 center(a + 0.9 * drand_r(), 0.2, b + 0.9 * drand_r());
      if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
        if (choose_mat < 0.8) { // diffuse
          //          list[i++] = new sphere(
          list[i++] = new moving_sphere(
              center, center + vec3(0, 0.5 * drand_r(), 0), 0.0, 1.0, 0.2,
              new lambertian(new constant_texture(
                  vec3(drand_r() * drand_r(), drand_r() * drand_r(),
                       drand_r() * drand_r()))));
        } else if (choose_mat < 0.95) { // metal
          list[i++] = new moving_sphere(
              center, center + vec3(-0.5f + drand_r(), -0.5f + drand_r(), 0.0),
              0.0, 1.0, 0.2,
              new metal(vec3(0.5 * (1 + drand_r()), 0.5 * (1 + drand_r()),
                             0.5 * (1 + drand_r())),
                        0.5 * drand_r()));
        } else { // glass
          list[i++] = new sphere(center, 0.2, new dielectric(1.5));
        }
      }
    }
  }

  list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
  list[i++] =
      new sphere(vec3(-4, 1, 0), 1.0,
                 new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
  list[i++] =
      new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

  return new bvh_node(list, i, 0, 1);
  //  return new hitable_list(list, i);
  //  return new hitable_list(list_bvh, k);
}

hitable *two_spheres() {
  texture *checker =
      new checker_texture(new constant_texture({0.1, 0.1, 0.1}), // white
                          new constant_texture({0.9, 0.9, 0.9})  // black
      );
  int n = 3;
  hitable **list = new hitable *[n];
  list[0] = new sphere(vec3(0, 10, 0), 10, new lambertian(checker));
  list[1] = new sphere(vec3(0, -10, 0), 10, new lambertian(checker));
  return new hitable_list(list, 2);
}

std::array<float, 256> perlin_noise::random_float_;
std::array<int, 256> perlin_noise::permute_x_;
std::array<int, 256> perlin_noise::permute_y_;
std::array<int, 256> perlin_noise::permute_z_;

float perlin_noise::noise(const vec3 &p) const {
  float u = p.x() - std::floor(p.x());
  float v = p.y() - std::floor(p.y());
  float w = p.z() - std::floor(p.z());

  auto valid_index = [](int &&index) {
    //    int tmp = index % 255;
    //    return tmp > 0 ? tmp : (255 + tmp);
    //  equals
    return index & 255;
  };
  int i = valid_index(std::floor(p.x()));
  int j = valid_index(std::floor(p.y()));
  int k = valid_index(std::floor(p.z()));

  float vertices[2][2][2]; // in a cube, tri-linear interpolate
  for (int ti = 0; ti < 2; ti++) {
    for (int tj = 0; tj < 2; tj++) {
      for (int tk = 0; tk < 2; tk++) {
        vertices[ti][tj][tk] = random_float_[permute_x_[valid_index(i + ti)] ^
                                             permute_y_[valid_index(j + tj)] ^
                                             permute_z_[valid_index(k + tk)]];
      }
    }
  }
  //  return trilinear_interpolate(vertices, u, v, w);
  return perlin_interpolate(vertices, u, v, w);

  // I guess it's a stable hash method to lookup float from random_float
  return random_float_[permute_x_[i] ^ permute_y_[j] ^ permute_z_[k]];
}

float perlin_noise::trilinear_interpolate(float vertex[2][2][2], float u,
                                          float v, float w) const {
  float accum = 0.0f;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        accum += (i * u + (1 - i) * (1 - u)) * (j * v + (1 - j) * (1 - v)) *
                 (k * w + (1 - k) * (1 - w)) * vertex[i][j][k];
      }
    }
  }
  return accum;
}

float perlin_noise::perlin_interpolate(float vertex[2][2][2], float u, float v,
                                       float w) const {
  auto trans = [](float &input) {
    return (6 * std::pow(input, 5) - 15 * std::pow(input, 4) +
            10 * std::pow(input, 3));
  };
  u = trans(u);
  v = trans(v);
  w = trans(w);
  float accum = 0.0f;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        accum += (i * u + (1 - i) * (1 - u)) * (j * v + (1 - j) * (1 - v)) *
                 (k * w + (1 - k) * (1 - w)) * vertex[i][j][k];
      }
    }
  }
  return accum;
}

hitable *two_perlin_spheres() {
  texture *perlin_texture = new perlin_noise_texture(10.0f);
  int n = 3;
  hitable **list = new hitable *[n];
  list[0] = new sphere(vec3(0, 2, 0), 2, new lambertian(perlin_texture));
  list[1] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(perlin_texture));
  return new hitable_list(list, 2);
}
