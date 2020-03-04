#include "camera.h"
#include "hitableList.h"
#include "ray.h"
#include "sphere.h"
#include "third_party/inipp.h"
#include "vec3.h"
#include <algorithm>
#include <condition_variable>
#include <deque>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <thread>

#ifdef __linux__
#include <cstdlib> //std::system
#endif

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
      float choose_mat = drand_r();
      vec3 center(a + 0.9 * drand_r(), 0.2, b + 0.9 * drand_r());
      if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
        if (choose_mat < 0.8) { // diffuse
          //          list[i++] = new sphere(
          list[i++] = new moving_sphere(
              center, center + vec3(0, 0.5 * drand_r(), 0), 0.0, 1.0, 0.2,
              new lambertian(vec3(drand_r() * drand_r(), drand_r() * drand_r(),
                                  drand_r() * drand_r())));
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
      new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
  list[i++] =
      new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

  hitable **list_bvh = new hitable *[10];
  int k = 0;
  return new bvh_node(list, i, 0, 1);
  //  return new hitable_list(list, i);
  //  return new hitable_list(list_bvh, k);
}
int main(int argc, char **argv) {
  std::string filename{"img.ppm"};
  std::ofstream ofs(filename, std::ios::out);

  std::ifstream config("config.ini");
  // default value
  int nx = 400;
  int ny = 200;
  int ns = 10; // anti-alisaing sample
  float aperture = 0.2; // the larger, the edge becomes more blurred
  float time0 = 0.0f;
  float time1 = 1.0f;
  // read value from config
  inipp::Ini<char> ini;
  ini.parse(config);
  inipp::extract(ini.sections["DEFAULT"]["width"], nx);
  inipp::extract(ini.sections["DEFAULT"]["height"], ny);
  inipp::extract(ini.sections["DEFAULT"]["sample"], ns);

  inipp::extract(ini.sections["BLUR"]["aperture"], aperture);

  inipp::extract(ini.sections["CAM_MOTION"]["start_time"], time0);
  inipp::extract(ini.sections["CAM_MOTION"]["end_time"], time1);

  ini.generate(std::cout);

  ofs << "P3\n" << nx << " " << ny << "\n255\n";

  hitable *world = random_scene();
  vec3 lookfrom = vec3(13, 2, 3);
  vec3 lookat = vec3(0, 0, 0);
  float dist_to_focus = (lookfrom - lookat).length();
  //  camera_with_blur cam(lookfrom, lookat, vec3(0, 1, 0), 90.0,
  //                       float(nx) / (float)ny, aperture, dist_to_focus);
  camera cam(lookfrom, lookat, vec3(0, 1, 0), 90.0, float(nx) / (float)ny,
             aperture, dist_to_focus, time0, time1);

  std::mutex mutex_;
  std::map<int, std::vector<int>> result;
  std::deque<std::thread> thread_vec;
  auto start = std::chrono::high_resolution_clock::now();
  for (int j = ny - 1; j >= 0; j--) {
    // naive thread pool
    //        if (thread_vec.size() >= 8) {
    //          std::for_each(thread_vec.begin(), thread_vec.end(),
    //                        [](std::thread &t) { t.join(); });
    //          thread_vec.clear();
    //        }
    thread_vec.emplace_back(
        [&](int index) {
          thread_local std::vector<int> row_colors;
          for (int i = 0; i < nx; i++) {
            vec3 col = vec3(0.0, 0.0, 0.0);
            for (int k = 0; k < ns; k++) {
              float u = ((float)i + drand_r()) / (float)nx;
              float v = ((float)index + drand_r()) / (float)ny;

              ray r = cam.get_ray(u, v);
              col += color(r, world, 0);
            }
            col /= float(ns);
            col = vec3(sqrt(col.r()), sqrt(col.g()), sqrt(col.b()));
            int red = int(255.99f * col.r());
            int green = int(255.99f * col.g());
            int blue = int(255.99f * col.b());
            row_colors.push_back(red);
            row_colors.push_back(green);
            row_colors.push_back(blue);
          }
          {
            std::lock_guard<std::mutex> lock(mutex_);
            result.insert({index, row_colors});
          }
        },
        j);
  }
  for (auto &i : thread_vec) {
    i.join();
  }
  for (int i = ny - 1; i >= 0; i--) {
    auto j = result.at(i);
    for (int k = 0; k < j.size(); k += 3) {
      ofs << j.at(k) << ' ' << j.at(k + 1) << ' ' << j.at(k + 2) << ' ' << '\n';
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << std::endl;
#ifdef __linux__
  // convert ppm to jpg
  // need imagemagic & Linux
  ofs.flush();
  std::string command = "convert " + filename + " " +
                        filename.substr(0, filename.find_first_of('.')) +
                        ".jpg";
  std::system(command.c_str());
#endif
  return 0;
}
