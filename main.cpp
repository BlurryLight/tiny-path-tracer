#include "camera.h"
#include "hitableList.h"
#include "ray.h"
#include "sphere.h"
#include "vec3.h"
#include <fstream>
#include <iostream>
#include <random>

vec3 color(const ray &r, hitable *world, int depth) {
  hit_record rec;
  if (world->hit(r, 0.001, MAXFLOAT, rec)) {
    ray scattered;
    vec3 attenuation;
    if (depth < 3 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
      return attenuation * color(scattered, world, depth + 1);
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
  list[0] =
      new sphere(vec3(0, -1000, 0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
  int i = 1;
  for (int a = -11; a < 11; a++) {
    for (int b = -11; b < 11; b++) {
      float choose_mat = drand48();
      vec3 center(a + 0.9 * drand48(), 0.2, b + 0.9 * drand48());
      if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
        if (choose_mat < 0.8) { // diffuse
          list[i++] = new sphere(
              center, 0.2,
              new lambertian(vec3(drand48() * drand48(), drand48() * drand48(),
                                  drand48() * drand48())));
        } else if (choose_mat < 0.95) { // metal
          list[i++] = new sphere(
              center, 0.2,
              new metal(vec3(0.5 * (1 + drand48()), 0.5 * (1 + drand48()),
                             0.5 * (1 + drand48())),
                        0.5 * drand48()));
        } else { // glass
          list[i++] = new sphere(center, 0.2, new dielectric(1.5));
        }
      }
    }
  }

  list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
  list[i++] =
      new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
  list[i++] =
      new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

  return new hitable_list(list, i);
}
int main(int argc, char **argv) {
  std::string filename{"img.ppm"};
  std::ofstream ofs(filename, std::ios::out);
  int nx = 400;
  int ny = 200;
  int ns = 100; // anti-alisaing sample
  ofs << "P3\n" << nx << " " << ny << "\n255\n";

  hitable *world = random_scene();
  vec3 lookfrom = vec3(9, 2, 3);
  vec3 lookat = vec3(0, 0, -1);
  float dist_to_focus = (lookfrom - lookat).length();
  float aperture = 0.1;
  camera_with_blur cam(lookfrom, lookat, vec3(0, 1, 0), 90.0,
                       float(nx) / (float)ny, aperture, dist_to_focus);
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      vec3 col = vec3(0.0, 0.0, 0.0);
      for (int k = 0; k < ns; k++) {
        float u = ((float)i + drand48()) / (float)nx;
        float v = ((float)j + drand48()) / (float)ny;

        ray r = cam.get_ray(u, v);
        col += color(r, world, 0);
      }
      col /= float(ns);
      col = vec3(sqrt(col.r()), sqrt(col.g()), sqrt(col.b()));
      int red = int(255.99f * col.r());
      int green = int(255.99f * col.g());
      int blue = int(255.99f * col.b());
      { ofs << red << ' ' << green << ' ' << blue << ' ' << '\n'; }
    }
  }
  return 0;
}
