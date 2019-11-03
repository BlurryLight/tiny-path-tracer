#include "camera.h"
#include "hitableList.h"
#include "ray.h"
#include "sphere.h"
#include "vec3.h"
#include <iostream>
#include <random>

vec3 color(const ray &r, hitable *world) {
  hit_record rec;
  if (world->hit(r, 0.0, MAXFLOAT, rec)) {
    return 0.5 *
           vec3(rec.normal.x() + 1, rec.normal.y() + 1, rec.normal.z() + 1);
  } else {

    vec3 unit_direction = unit_vector(r.direction());
    float t = (unit_direction.y() + 1.0) * 0.5; // clamp (-1,1) to (0,1)
    return (1 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
  }
  // linear interperation
  // blend white with blue
}

int main()
{
  int nx = 200;
  int ny = 100;
  int ns = 100;
  std::cout << "P3\n" << nx << " " << ny << "\n255\n";
  hitable *list[2];
  list[0] = new sphere(vec3(0, 0, -1), 0.5);
  list[1] = new sphere(vec3(0, -100.5, -1), 100);

  hitable *world = new hitable_list(list, 2);
  camera cam;
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {

      vec3 col = vec3(0.0, 0.0, 0.0);
      for (int k = 0; k < ns; k++) {

        float u = ((float)i + drand48()) / (float)nx;
        float v = ((float)j + drand48()) / (float)ny;

        ray r = cam.get_ray(u, v);
        col += color(r, world);
      }
      col /= float(ns);

      int red = int(255.99f * col.r());
      int green = int(255.99f * col.g());
      int blue = int(255.99f * col.b());
      std::cout << red << ' ' << green << ' ' << blue << ' ' << '\n';
    }
  }
  return 0;
}
