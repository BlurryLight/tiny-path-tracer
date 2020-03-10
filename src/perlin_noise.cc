#include "perlin_noise.h"
#include "utils.h"
perlin_noise::perlin_noise() {
  std::generate(random_vec3_.begin(), random_vec3_.end(), []() -> vec3 {
    float f1 = static_cast<float>(2 * drand_r() - 1);
    float f2 = static_cast<float>(2 * drand_r() - 1);
    float f3 = static_cast<float>(2 * drand_r() - 1);
    return unit_vector(vec3(f1, f2, f3));
  });
  for (int i = 0; i != permute_x_.size(); i++) {
    permute_x_[i] = i;
    permute_y_[i] = i;
    permute_z_[i] = i;
  }
  auto seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  auto engine = std::default_random_engine(seed);
  std::shuffle(permute_x_.begin(), permute_x_.end(), engine);
  std::shuffle(permute_y_.begin(), permute_y_.end(), engine);
  std::shuffle(permute_z_.begin(), permute_z_.end(), engine);
}
