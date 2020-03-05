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

  return vec3(1, 1, 1) * noise_.noise(scale_ * p);
}
