#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
struct vec3 {
  float e[3];
  vec3() {}
  vec3(float e0, float e1, float e2) {
    e[0] = e0;
    e[1] = e1;
    e[2] = e2;
  }
  inline float x() { return e[0]; }
  inline float y() { return e[1]; }
  inline float z() { return e[2]; }

  inline float r() { return e[0]; }
  inline float g() { return e[1]; }
  inline float b() { return e[2]; }

  inline const vec3 &operator+() const { return *this; }
  inline vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }

  inline float operator[](int index) const { return e[index]; }
  inline float &operator[](int index) { return e[index]; }

  inline vec3 &operator+=(const vec3 &v2);
  inline vec3 &operator-=(const vec3 &v2);
  inline vec3 &operator*=(const vec3 &v2);
  inline vec3 &operator/=(const vec3 &v2);
  inline vec3 &operator*=(const float t);
  inline vec3 &operator/=(const float t);

  inline float length() const {
    return std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  }
  inline float squared_length() const {
    return (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  }
  inline void make_unit_vector();
  inline void normalize() { make_unit_vector(); }
};

inline vec3 operator+(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
}
inline vec3 operator-(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
}
inline vec3 operator*(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
}
inline vec3 operator/(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
}

inline vec3 operator*(const vec3 &v1, float t) {
  return vec3(v1.e[0] * t, v1.e[1] * t, v1.e[2] * t);
}
inline vec3 operator/(const vec3 &v1, float t) {
  return vec3(v1.e[0] / t, v1.e[1] / t, v1.e[2] / t);
}
inline vec3 operator*(float t, const vec3 &v1) {
  return vec3(v1.e[0] * t, v1.e[1] * t, v1.e[2] * t);
}

inline float dot(const vec3 &v1, const vec3 &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}
inline vec3 cross(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[1] * v2.e[2] - v2.e[1] * v1.e[2],
              v1.e[2] * v2.e[0] - v2.e[2] * v1.e[0],
              v1.e[0] * v2.e[1] - v2.e[0] * v1.e[1]);
}

inline std::istream &operator>>(std::istream &is, vec3 &vec) {
  is >> vec.e[0] >> vec.e[1] >> vec.e[2];
  return is;
}

inline std::ostream &operator<<(std::ostream &os, const vec3 &vec) {
  os << vec.e[0] << " " << vec.e[1] << " " << vec.e[2];
  return os;
}
inline void vec3::make_unit_vector() {
  float k = 1.0 / this->length();
  e[0] *= k;
  e[1] *= k;
  e[2] *= k;
}

vec3 &vec3::operator+=(const vec3 &v2) {
  e[0] += v2.e[0];
  e[1] += v2.e[1];
  e[2] += v2.e[2];
  return *this;
}

vec3 &vec3::operator-=(const vec3 &v2) {
  e[0] -= v2.e[0];
  e[1] -= v2.e[1];
  e[2] -= v2.e[2];
  return *this;
}

vec3 &vec3::operator*=(const vec3 &v2) {
  e[0] *= v2.e[0];
  e[1] *= v2.e[1];
  e[2] *= v2.e[2];
  return *this;
}

vec3 &vec3::operator/=(const vec3 &v2) {
  e[0] /= v2.e[0];
  e[1] /= v2.e[1];
  e[2] /= v2.e[2];
  return *this;
}

vec3 &vec3::operator/=(const float t) {
  e[0] /= t;
  e[1] /= t;
  e[2] /= t;
  return *this;
}
vec3 &vec3::operator*=(const float t) {
  e[0] *= t;
  e[1] *= t;
  e[2] *= t;
  return *this;
}

inline vec3 normalize(vec3 v) { return v / v.length(); }
inline vec3 unit_vector(vec3 v) { return normalize(v); }

#endif // VEC3_H
