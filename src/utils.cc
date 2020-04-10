#include "utils.h"
#include "hitable.h"
#include "hitable_list.h"
#include "material.h"
#include "perlin_noise.h"
#include "ray.h"
#include "rect_box.h"
#include <memory>
#include <random>
#include <sphere.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
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
  } while (p.length() >= 1.0);
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
vec3 color_shadow(const ray &r, hitable *world, hitable *light_shape, int depth,
           int max_depth) {
  hit_record hit_rec;
  if (world->hit(r, 0.001, std::numeric_limits<float>::max(), hit_rec)) {
    scatter_record scatter_rec;
    vec3 emitted = hit_rec.mat_ptr->emitted(r, hit_rec, hit_rec.u, hit_rec.v,
                                            hit_rec.point);
    if (depth < max_depth &&
        hit_rec.mat_ptr->scatter(r, hit_rec, scatter_rec)) {
      if (scatter_rec.is_specular) {
        //金属直接反射，不用考虑pdf
        return scatter_rec.attenuation * color(scatter_rec.specular_ray, world,
                                               light_shape, depth + 1,
                                               max_depth);
      }
      vec3 Light_color{20,20,20};
      vec3 L_dir{0,0,0};
      auto shadow_ray = ray(hit_rec.point,light_shape->random(hit_rec.point));
      hit_record shadow_rec;
      if (world->hit(shadow_ray, 0.001, std::numeric_limits<float>::max(), shadow_rec)) {
        if(std::abs(shadow_rec.point.y() - 298.0f) < 1.0f)
            {
          L_dir = Light_color * dot(normalize(shadow_ray.direction()),hit_rec.normal) / light_shape->pdf_value(hit_rec.point,shadow_ray.direction());
          if (std::isnan(L_dir.x()) || std::isnan(L_dir.y()) ||
              std::isnan(L_dir.z()))
          {
            L_dir = vec3(0.0f,0.0f,0.0f);
          }
        }

      }

      auto scattered = ray(hit_rec.point,scatter_rec.pdf_ptr->generate());
          return emitted + scatter_rec.attenuation * (L_dir +
                        color(scattered,world,light_shape,depth+1,max_depth));
    }
    else {
      return emitted;
    }
  } else {
    return vec3(0, 0, 0);
    //    vec3 unit_direction = unit_vector(r.direction());
    //    float t = (unit_direction.y() + 1.0) * 0.5; // clamp (-1,1) to (0,1)
    //    return (0.1) * ((1 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5,
    //    0.7, 1.0));
  }
  // linear interperation
  // blend white with blue

}

vec3 color(const ray &r, hitable *world, hitable *light_shape, int depth,
           int max_depth) {
  hit_record hit_rec;
  if (world->hit(r, 0.001, std::numeric_limits<float>::max(), hit_rec)) {
    scatter_record scatter_rec;
    vec3 emitted = hit_rec.mat_ptr->emitted(r, hit_rec, hit_rec.u, hit_rec.v,
                                            hit_rec.point);
    if (depth < max_depth &&
        hit_rec.mat_ptr->scatter(r, hit_rec, scatter_rec)) {
      if (scatter_rec.is_specular) {
        //金属直接反射，不用考虑pdf
        return scatter_rec.attenuation * color(scatter_rec.specular_ray, world,
                                               light_shape, depth + 1,
                                               max_depth);
      }
      hitable_pdf p0(light_shape, hit_rec.point);
      mixture_pdf p(&p0, scatter_rec.pdf_ptr.get());
      ray scattered = ray(hit_rec.point, p.generate(), r.time());
      float pdf_value = p.value(scattered.direction());
      return emitted +
             scatter_rec.attenuation *
                 hit_rec.mat_ptr->scattering_pdf(r, hit_rec, scattered) *
                 color(scattered, world, light_shape, depth + 1, max_depth) /
                 pdf_value;
    } else {
      return emitted;
    }
  } else {
    return vec3(0, 0, 0);
    //    vec3 unit_direction = unit_vector(r.direction());
    //    float t = (unit_direction.y() + 1.0) * 0.5; // clamp (-1,1) to (0,1)
    //    return (0.1) * ((1 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5,
    //    0.7, 1.0));
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

hitable *two_checker_spheres() {
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

std::array<vec3, 256> perlin_noise::random_vec3_;
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

  vec3 vertices_vec[2][2][2];
  for (int ti = 0; ti < 2; ti++) {
    for (int tj = 0; tj < 2; tj++) {
      for (int tk = 0; tk < 2; tk++) {
        vec3 tmp = random_vec3_[permute_x_[valid_index(i + ti)] ^
                                permute_y_[valid_index(j + tj)] ^
                                permute_z_[valid_index(k + tk)]];
        vertices_vec[ti][tj][tk] = tmp;
      }
    }
  }
  //  return trilinear_interpolate(vertices, u, v, w);
  return perlin_interpolate(vertices_vec, u, v, w);

  // I guess it's a stable hash method to lookup float from random_float
  //  return random_float_[permute_x_[i] ^ permute_y_[j] ^ permute_z_[k]];
}

float perlin_noise::turb(const vec3 &p, int depth) const {
  float accum = 0;
  vec3 tmp = p;
  float weight = 1.0;
  for (int i = 0; i < depth; i++) {
    accum += weight * noise(tmp);
    weight *= 0.5;
    tmp *= 2;
  }
  return std::abs(accum);
}

float perlin_noise::perlin_interpolate(vec3 vertex[2][2][2], float u, float v,
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
        vec3 weight_v(u - i, v - j, w - k);
        accum += (i * u + (1 - i) * (1 - u)) * (j * v + (1 - j) * (1 - v)) *
                 (k * w + (1 - k) * (1 - w)) * dot(vertex[i][j][k], weight_v);
      }
    }
  }
  return std::abs(accum);
}

hitable *two_perlin_spheres() {
  texture *perlin_texture = new perlin_noise_texture(2.0f);
  int n = 3;
  hitable **list = new hitable *[n];
  list[0] = new sphere(vec3(0, 2, 0), 2, new lambertian(perlin_texture));
  list[1] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(perlin_texture));
  return new hitable_list(list, 2);
}

unsigned char *load_image_texture(std::string filename, int &width, int &height,
                                  int &channels) {
  auto data = stbi_load(filename.c_str(), &width, &height, &channels, 0);
  return data;
}

hitable *light_spheres() {
  texture *light_texture = (new constant_texture(vec3(4, 4, 4)));
  int n = 10;
  hitable **list = new hitable *[n];
  list[0] =
      new sphere(vec3(-1, 1, 0), 1,
                 new lambertian(new constant_texture(vec3(0.3, 0.4, 0.5))));
  list[1] = new sphere(vec3(-3, 1, 2), 1, new diffuse_light(light_texture));
  list[2] = new sphere(vec3(0, -1000, 0), 1000,
                       new lambertian(new perlin_noise_texture(4.0f)));
  list[3] = new sphere(vec3(-3, 1, -2), 1, new diffuse_light(light_texture));
  list[4] = new xy_rect(3, 5, 1, 3, -2, new diffuse_light(light_texture));
  return new hitable_list(list, 5);
}

hitable *sphere_cornell_box() {
  hitable **list = new hitable *[100];
  int i = 0;
  material *red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
  material *white =
      new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
  material *green =
      new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
  material *light = new diffuse_light(new constant_texture(vec3(20, 20, 20)));

  list[i++] = new sphere(vec3(-1e5, 0, 0), 1e5 - 300, green);
  list[i++] = new sphere(vec3(1e5, 0, 0), 1e5 - 300, red);
  list[i++] = new sphere(vec3(0, 1e5, 0), 1e5 - 300, white);
  list[i++] = new sphere(vec3(0, -1e5, 0), 1e5 - 300, white);
  list[i++] = new sphere(vec3(0, 0, -1e5), 1e5 - 300, white);
  list[i++] = new xz_rect(-150, 150, -150, 150, 300, light);
  list[i++] = new sphere(vec3(-150, -200, -100), 100, new dielectric(1.5));
  list[i++] = new sphere(vec3(150, -200, 100), 100,
                         new metal(vec3(0.7, 0.6, 0.5), 0.0));
  //  list[i++] = new flip_normal(new yz_rect(0, 555, 0, 555, 555, green));
  //  list[i++] = new yz_rect(0, 555, 0, 555, 0, red);
  //  list[i++] = new xz_rect(213, 343, 227, 332, 554, light);
  //  list[i++] = new flip_normal(new xz_rect(0, 555, 0, 555, 555, white));
  //  list[i++] = new xz_rect(0, 555, 0, 555, 0, white);
  //  list[i++] = new flip_normal(new xy_rect(0, 555, 0, 555, 555, white));

  //  return new hitable_list(list, i);
  return new bvh_node(list, i, 0, 0);
}

hitable *cornell_box() {
  hitable **list = new hitable *[100];
  int i = 0;
  material *red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
  material *white =
      new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
  material *green =
      new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
  material *light = new diffuse_light(new constant_texture(vec3(20, 20, 20)));

  list[i++] = new yz_rect(-300, 300, -300, 300, -300, green); // left
  list[i++] =
      new flip_normal(new yz_rect(-300, 300, -300, 300, 300, red));    // right
  list[i++] = new xz_rect(-300, 300, -300, 300, -300, white);          // bottom
  list[i++] =
      new flip_normal(new xz_rect(-300, 300, -300, 300, 300, white)); // top
  list[i++] = new xy_rect(-300, 300, -300, 300, -300, white);
  list[i++] = new flip_normal(new xz_rect(-100, 100, -150, -50, 298, light));

  vec3 rect_box_corner = vec3(-200, -300, -100);
  material *aluminum = new metal(vec3(0.8, 0.85, 0.88), 0.0);
  list[i++] = new translate(
      new rotate_y(new box(vec3(0, 0, 0), vec3(200, 350, 75), aluminum), 35.0f),
      rect_box_corner); // rectangle
  vec3 square_box_corner = vec3(30, -300, -50);
  list[i++] = new translate(
      new rotate_y(new box(vec3(0, 0, 0), vec3(180, 180, 180), white), -25.0f),
      square_box_corner); // square
  list[i++] = new sphere(vec3(120, -50, 40), 70, new dielectric(1.5));

  //  return new hitable_list(list, i);
  return new bvh_node(list, i, 0, 0);
}

hitable *cornell_box_smoke() {
  hitable **list = new hitable *[100];
  int i = 0;
  material *red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
  material *white =
      new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
  material *green =
      new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
  material *light = new diffuse_light(new constant_texture(vec3(20, 20, 20)));

  list[i++] = new yz_rect(-300, 300, -300, 300, -300, green); // left
  list[i++] =
      new flip_normal(new yz_rect(-300, 300, -300, 300, 300, red));    // right
  list[i++] = new xz_rect(-300, 300, -300, 300, -300, white);          // bottom
  list[i++] =
      new flip_normal(new xz_rect(-300, 300, -300, 300, 300, white)); // top
  list[i++] = new xy_rect(-300, 300, -300, 300, -300, white);
  list[i++] = new flip_normal(new xz_rect(-100, 100, -150, -50, 298, light));

  vec3 rect_box_corner = vec3(-200, -300, -100);
  auto rec_box = new translate(
      new rotate_y(new box(vec3(0, 0, 0), vec3(200, 350, 75), white), 45.0f),
      rect_box_corner); // rectangle
  vec3 square_box_corner = vec3(30, -300, -50);
  auto square_box = new translate(
      new rotate_y(new box(vec3(0, 0, 0), vec3(180, 180, 180), white), -15.0f),
      square_box_corner); // square
  list[i++] = new constant_medium(rec_box, 0.05,
                                  new constant_texture(vec3(1.0, 1.0, 1.0)));
  list[i++] = new constant_medium(square_box, 0.01,
                                  new constant_texture(vec3(0.1, 0.0, 0.0)));
  hitable *mist = new sphere(vec3(0, 0, 0), 1000, new dielectric(1.5));
  list[i++] = new constant_medium(mist, 0.0001,
                                  new constant_texture(vec3(1.0, 1.0, 1.0)));

  return new hitable_list(list, i);
}

hitable *oneweek_final() {
  int nb = 20;
  hitable **list = new hitable *[30];
  hitable **boxlist = new hitable *[10000];
  hitable **boxlist2 = new hitable *[10000];
  material *white =
      new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
  material *ground =
      new lambertian(new constant_texture(vec3(0.48, 0.83, 0.53)));
  int b = 0;
  for (int i = 0; i < nb; i++) {
    for (int j = 0; j < nb; j++) {
      float w = 100;
      float x0 = -1000 + i * w;
      float z0 = -1000 + j * w;
      float y0 = 0;
      float x1 = x0 + w;
      float y1 = 100 * (drand_r() + 0.01);
      float z1 = z0 + w;
      boxlist[b++] = new box(vec3(x0, y0, z0), vec3(x1, y1, z1), ground);
    }
  }
  int l = 0;
  list[l++] = new bvh_node(boxlist, b, 0, 1);
  material *light = new diffuse_light(new constant_texture(vec3(7, 7, 7)));
  list[l++] = new xz_rect(123, 423, 147, 412, 554, light);
  vec3 center(400, 400, 200);
  list[l++] = new moving_sphere(
      center, center + vec3(30, 0, 0), 0, 1, 50,
      new lambertian(new constant_texture(vec3(0.7, 0.3, 0.1))));
  list[l++] = new sphere(vec3(260, 150, 45), 50, new dielectric(1.5));
  list[l++] =
      new sphere(vec3(0, 150, 145), 50, new metal(vec3(0.8, 0.8, 0.9), 10.0));
  hitable *boundary = new sphere(vec3(360, 150, 145), 70, new dielectric(1.5));
  list[l++] = boundary;
  list[l++] = new constant_medium(boundary, 0.2,
                                  new constant_texture(vec3(0.2, 0.4, 0.9)));
  boundary = new sphere(vec3(0, 0, 0), 5000, new dielectric(1.5));
  list[l++] = new constant_medium(boundary, 0.0001,
                                  new constant_texture(vec3(1.0, 1.0, 1.0)));
  int nx, ny, nn;
  unsigned char *tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);
  material *emat = new lambertian(new image_texture(tex_data, nx, ny));
  list[l++] = new sphere(vec3(400, 200, 400), 100, emat);
  texture *pertext = new perlin_noise_texture(0.1);
  list[l++] = new sphere(vec3(220, 280, 300), 80, new lambertian(pertext));
  int ns = 1000;
  for (int j = 0; j < ns; j++) {
    boxlist2[j] = new sphere(
        vec3(165 * drand_r(), 165 * drand_r(), 165 * drand_r()), 10, white);
  }
  list[l++] =
      new translate(new rotate_y(new bvh_node(boxlist2, ns, 0.0, 1.0), 15),
                    vec3(-100, 270, 395));
  return new hitable_list(list, l);
}

vec3 random_on_sphere() {
  // some magic here
  // keyword: uniform sample on a sphere
  float r1 = drand_r();
  float r2 = drand_r();
  float x = std::cos(2 * M_PI * r1) * 2 * std::sqrt(r2 * (1 - r2));
  float y = std::sin(2 * M_PI * r1) * 2 * std::sqrt(r2 * (1 - r2));
  float z = 1 - 2 * r2;
  return vec3(x, y, z);
}

vec3 random_on_hemisphere() {
  float r1 = drand_r();
  float r2 = drand_r();
  float phi = 2 * M_PI * r1;
  float x = std::cos(phi) * std::sqrt(r2);
  float y = std::sin(phi) * std::sqrt(r2);
  float z = std::sqrt(1 - r2);
  return vec3(x, y, z);
}

void onb::build_from_w(const vec3 &normal) {
  axis_[2] = normal;
  vec3 tmp;
  if (std::abs(normal.x()) > 0.9) // is normal x-axis?
  {
    // normal is x-axis ,so let tmp be y-axis
    tmp = vec3(0, 1, 0);
  } else {
    tmp = vec3(1, 0, 0);
  }

  axis_[1] = unit_vector(cross(normal, tmp));
  axis_[0] = cross(v(), w());
}

float hitable_pdf::value(const vec3 &direction) const {
  return ptr_->pdf_value(origin_, direction);
}

vec3 hitable_pdf::generate() const { return ptr_->random(origin_); }
