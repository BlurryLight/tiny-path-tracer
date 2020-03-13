#ifndef HISTABLELIST_H
#define HISTABLELIST_H
#include "hitable.h"
class hitable_list : public hitable {
public:
  hitable **list_;
  int list_size_;
  hitable_list() {}
  hitable_list(hitable **l, int n) {
    list_ = l;
    list_size_ = n;
  }
  virtual bool hit(const ray &r, float t_min, float t_max,
                   hit_record &rec) const override;
  virtual bool bounding_box(float t0, float t1, AABB &box) const override;
  virtual float pdf_value(const vec3 &origin,
                          const vec3 &direction) const override;
  virtual vec3 random(const vec3 &origin) const override;
};

#endif // HISTABLELIST_H
