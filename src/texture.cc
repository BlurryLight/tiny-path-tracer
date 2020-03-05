#include "texture.h"
#include "utils.h"

vec3 checker_texture::value(float u, float v, const vec3 &p) const {
  //(0,0,0) = positive -> even
  //(0.2,0.2,0.2) = sin(2) * sin(2) * sin(2) = positive ->even
  //(0.4,0.4,0.4) = sin(4) * sin(4) * sin(4)  = negative -> odd texture
  float sin_sign =
      std::sin(10 * p.x()) * std::sin(10 * p.y()) * std::sin(10 * p.z());
  if (std::isless(sin_sign, 0.0f)) // negative
  {
    return odd_->value(u, v, p);
  } else {
    return even_->value(u, v, p);
  }
}

vec3 perlin_noise_texture::value(float u, float v, const vec3 &p) const {

  //  return vec3(1, 1, 1) * noise_.noise(scale_ * p);

  // with turb
  return vec3(1, 1, 1) * 0.5 *
         (1 + std::sin(scale_ * p.z() + 10 * noise_.turb(p)));
}

vec3 image_texture::value(float u, float v, const vec3 &p) const {
  int i = u * width_;
  int j = (1 - v) * height_;
  if (i < 0)
    i = 0;
  if (i > width_ - 1)
    i = width_ - 1;
  if (j > height_ - 1)
    j = height_ - 1;
  if (j < 0)
    j = 0;
  float r = int(data_[3 * i + 3 * width_ * j]) / 255.0f;
  float g = int(data_[3 * i + 3 * width_ * j + 1]) / 255.0f;
  float b = int(data_[3 * i + 3 * width_ * j + 2]) / 255.0f;
  return vec3(r, g, b);
}
