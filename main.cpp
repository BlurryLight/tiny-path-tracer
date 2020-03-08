#include "camera.h"
#include "ray.h"
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
#include <sphere.h>
#include <texture.h>
#include <thread>
#include <unordered_map>

#ifdef __linux__
#include <cstdlib> //std::system
#endif

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
  int sample_max_recurse_depth = 50;
  float fov = 90.0f;
  int bonus_pic = 10;
  // read value from config
  inipp::Ini<char> ini;
  ini.parse(config);
  inipp::extract(ini.sections["DEFAULT"]["width"], nx);
  inipp::extract(ini.sections["DEFAULT"]["height"], ny);
  inipp::extract(ini.sections["DEFAULT"]["sample"], ns);
  inipp::extract(ini.sections["DEFAULT"]["recur_depth"],
                 sample_max_recurse_depth);
  inipp::extract(ini.sections["DEFAULT"]["fov"], fov);
  inipp::extract(ini.sections["DEFAULT"]["fov"], bonus_pic);

  inipp::extract(ini.sections["BLUR"]["aperture"], aperture);

  inipp::extract(ini.sections["CAM_MOTION"]["start_time"], time0);
  inipp::extract(ini.sections["CAM_MOTION"]["end_time"], time1);

  ini.generate(std::cout);

  std::cout << "<===========>" << std::endl;
#ifdef DEBUG_MODE
  std::cout << "IN DEBUG MODE" << std::endl;
#else
  std::cout << "IN RELEASE MODE" << std::endl;
#endif
  std::cout << "<===========>" << std::endl;
  ofs << "P3\n" << nx << " " << ny << "\n255\n";

  //  hitable *world = sphere_cornell_box();
  hitable *world = cornell_box();
  //  hitable *world = random_scene();
  //  hitable *world = light_spheres();
  //  hitable *world = two_perlin_spheres();
  //  int width, height, channels;
  //  auto data = load_image_texture("earthmap.jpg", width, height, channels);
  //  hitable *world = new sphere(
  //      {0, 0, 0}, 3, new lambertian(new image_texture(data, width, height)));
  //  vec3 lookfrom = vec3(13, 2, 3);
  //  vec3 lookat = vec3(0, 0, 0);
  //  float dist_to_focus = (lookfrom - lookat).length();

  // cornel box
  vec3 lookfrom = vec3(0, 0, 800);
  vec3 lookat = vec3(0, 0, 0);
  float dist_to_focus = 10.0f;
  camera cam(lookfrom, lookat, vec3(0, 1, 0), fov, float(nx) / (float)ny,
             aperture, dist_to_focus, time0, time1);

  std::mutex mutex_;
  std::mutex mutex_2;
  std::map<int, std::vector<int>> result;
  std::unordered_map<std::string, std::vector<vec3>> pixel_sample_cols;
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
            std::vector<vec3> sample_cols;
            vec3 col = vec3(0.0, 0.0, 0.0);
            for (int k = 0; k < ns; k++) {
              float u = ((float)i + drand_r()) / (float)nx;
              float v = ((float)index + drand_r()) / (float)ny;

              ray r = cam.get_ray(u, v);
              auto tmp = color(r, world, 0, sample_max_recurse_depth);
              col += tmp;
              sample_cols.push_back(tmp);
            }
            col /= float(ns);
            col = vec3(sqrt(col.r()), sqrt(col.g()), sqrt(col.b()));
            int red = int(255.99f * col.r());
            int green = int(255.99f * col.g());
            int blue = int(255.99f * col.b());
            row_colors.push_back(red);
            row_colors.push_back(green);
            row_colors.push_back(blue);

            {
              std::lock_guard<std::mutex> lock(mutex_2);
              auto key = std::to_string(i) + "+" + std::to_string(index);
              pixel_sample_cols.insert({key, sample_cols});
            }
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
  auto valid_rgb = [](int rgb) {
    if (rgb < 0)
      rgb = 0;
    if (rgb > 255)
      rgb = 255;
    return rgb;
  };
  for (int i = ny - 1; i >= 0; i--) {
    auto j = result.at(i);
    for (int k = 0; k < j.size(); k += 3) {
      ofs << valid_rgb(j.at(k)) << ' ' << valid_rgb(j.at(k + 1)) << ' '
          << valid_rgb(j.at(k + 2)) << ' ' << '\n';
    }
  }
  std::cout << pixel_sample_cols.find("0+0")->second.size();
  std::cout.flush();

  int sample_vec_slice = ns / bonus_pic;
  for (int k = 0; k < bonus_pic; k++) {
    std::ofstream sample_ofs{"img_" + std::to_string(k) + ".ppm"};
    sample_ofs << "P3\n" << nx << " " << ny << "\n255\n";
    for (int j = ny - 1; j >= 0; j--) {
      for (int i = 0; i < nx; i++) {

        auto key = std::to_string(i) + "+" + std::to_string(j);
        std::vector<vec3> sample_vec = pixel_sample_cols.find(key)->second;
        std::vector<vec3> target_vec{sample_vec.cbegin(),
                                     sample_vec.cbegin() +
                                         sample_vec_slice * (k + 1)};

        vec3 col{0, 0, 0};
        for (const auto &it : target_vec) {
          col += it;
        }
        col /= float(sample_vec_slice * (k + 1));
        col = vec3(sqrt(col.r()), sqrt(col.g()), sqrt(col.b()));
        int red = int(255.99f * col.r());
        int green = int(255.99f * col.g());
        int blue = int(255.99f * col.b());
        sample_ofs << valid_rgb(red) << ' ' << valid_rgb(green) << ' '
                   << valid_rgb(blue) << ' ';
      }
    }
  }
  auto end = std::chrono::high_resolution_clock::now();

  // clang-format off
  std::cout << "time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() /1000.0f
            << " s" << std::endl;
#ifdef __linux__
  // convert ppm to jpg
  // need imagemagick & Linux
  ofs.flush();
  std::string command = "convert " + filename + " " +
                        filename.substr(0, filename.find_first_of('.')) +
                        ".jpg";
  std::system(command.c_str());
#endif
  return 0;
}
