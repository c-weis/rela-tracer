// Copyright 2023 Christoph Weis
// Basic 3- and 4-vector maths library
// Crucially includes Lorentz transforms
#include "include/math.h"

#include <cmath>
#include <iostream>
#include <random>

/*
    3-vector lib
*/

std::ostream &operator<<(std::ostream &out, const Vec3 &vec) {
  return out << "(" << vec.x << " " << vec.y << " " << vec.z << ")";
}

Vec3 operator+(const Vec3 &left, const Vec3 &right) {
  return Vec3(left.x + right.x, left.y + right.y, left.z + right.z);
}

Vec3 operator-(const Vec3 &left, const Vec3 &right) {
  return Vec3(left.x - right.x, left.y - right.y, left.z - right.z);
}

Vec3 operator-(const Vec3 &vec) { return Vec3(-vec.x, -vec.y, -vec.z); }

Vec3 operator*(const Vec3 &vec, float scalar) {
  return Vec3(vec.x * scalar, vec.y * scalar, vec.z * scalar);
}

Vec3 operator*(float scalar, const Vec3 &vec) { return vec * scalar; }

Vec3 operator/(const Vec3 &vec, float scalar) {
  return Vec3(vec.x / scalar, vec.y / scalar, vec.z / scalar);
}

Vec3 Vec3::NormalizedNonzero() const { return *this / Norm(); }

OptionalVec3 Vec3::Normalized() const {
  float norm = this->Norm();
  if (norm < kDivisionEpsilon) {
    return {};
  }
  return *this / norm;
}

float Vec3::NormSq() const { return Dot3(*this); }

float Vec3::Norm() const { return sqrtf(NormSq()); }

float Vec3::Dot3(const Vec3 &other) const {
  return x * other.x + y * other.y + z * other.z;
}

float Dot3(const Vec3 &left, const Vec3 &right) { return left.Dot3(right); }

Vec3 operator^(const Vec3 &left, const Vec3 &right) {
  return Cross(left, right);
}

float Vec3::CosAngleTo(const Vec3 &other) const {
  return Dot3(*this) / (Norm() * other.Norm());
}

Vec3 Vec3::Cross(const Vec3 &other) const {
  return Vec3(y * other.z - z * other.y, z * other.x - x * other.z,
              x * other.y - y * other.x);
}

Vec3 Cross(const Vec3 &left, const Vec3 &right) { return left.Cross(right); }

Vec3 NormalTo(const Vec3 &left, const Vec3 &right) {
  return Cross(left, right).Normalized().value_or(kDefaultNormal);
}

void Vec3::NormalizeNonzero() {
  float _norm = Norm();
  x /= _norm;
  y /= _norm;
  z /= _norm;
}

bool Vec3::Normalize() {
  /*
      Returns a bool indicating whether normalisation
      succeeded.
  */
  float _norm = Norm();
  if (_norm > kDivisionEpsilon) {
    x /= _norm;
    y /= _norm;
    z /= _norm;
    return true;
  }
  return false;
}

float Vec3::Gamma() const {
  // Return gamma for this considered
  // as a (relative) velocity
  return 1 / GammaInv();
}

float Vec3::GammaInv() const {
  // Return 1/gamma for this considered
  // as a (relative) velocity
  return sqrtf(1 - NormSq());
}

// Transform velocity. 
// Rf velocity is assumed to be sub-lightspeed(!)
Vec3 Vec3::VelTransformedToFrame(const Vec3 &vel) const {
  float gamma = vel.Gamma();
  float prod = Dot3(vel);
  return 1 / (1 - prod) *
         (*this / gamma - vel + gamma / (gamma + 1) * prod * vel);
}

Vec3 Vec3::VelTransformedFromFrame(const Vec3 &vel) const {
  return VelTransformedToFrame(-vel);
}

// Reflects vector in plane given a normal vector
// Normal is presumed to be normalised!
Vec3 Vec3::Reflect(const Vec3 &normal) const {
  return *this - 2 * Dot3(normal) * normal;
}

/*
    4-vector lib
*/

std::ostream &operator<<(std::ostream &out, const Vec4 &vec) {
  return out << "(" << vec.t << " " << vec.r << ")";
}

Vec4 operator+(const Vec4 &left, const Vec4 &right) {
  return Vec4(left.t + right.t, left.r + right.r);
}

Vec4 operator-(const Vec4 &left, const Vec4 &right) {
  return Vec4(left.t - right.t, left.r - right.r);
}

Vec4 operator-(const Vec4 &vec) { return Vec4(-vec.t, -vec.r); }

Vec4 operator*(const Vec4 &vec, float scalar) {
  return Vec4(vec.t * scalar, vec.r * scalar);
}

Vec4 operator*(float scalar, const Vec4 &vec) { return vec * scalar; }

Vec4 operator/(const Vec4 &vec, float scalar) {
  return Vec4(vec.t / scalar, vec.r / scalar);
}

float Vec4::Dot4(const Vec4 &other) const {
  // Minkowski metric
  return -t * other.t + r.Dot3(other.r);
}

float Dot4(const Vec4 &left, const Vec4 &right) { return left.Dot4(right); }

float Vec4::MinkowskiNormSq() const { return Dot4(*this); }

float Vec4::SpaceNormSq() const { return r.NormSq(); }

float Vec4::MinkowskiNorm() const { return sqrtf(MinkowskiNormSq()); }

float Vec4::SpaceNorm() const { return sqrtf(SpaceNormSq()); }

Vec4 Vec4::TransformedToFrame(const Vec3 &vel) const {
  // Lorentz boost position 4-vector
  float gamma = vel.Gamma();

  OptionalVec3 o_n = vel.Normalized();
  if (!o_n.has_value()) return Vec4(*this);
  Vec3 n = o_n.value();

  return Vec4(gamma * (t - Dot3(vel, r)),
              r + (gamma - 1) * Dot3(r, n) * n - gamma * t * vel);
}

Vec4 Vec4::TransformedFromFrame(const Vec3 &vel) const {
  // Lorentz de-boost position 4-vector
  return TransformedToFrame(-vel);
}

Vec4 Vec4::TransformedToFrame(const Line &other) const {
  return (*this - other.origin).TransformedToFrame(other.vel);
}

Vec4 Vec4::TransformedFromFrame(const Line &other) const {
  return TransformedFromFrame(other.vel) + other.origin;
}

/*
    Line through spacetime
*/
std::ostream &operator<<(std::ostream &out, const Line &line) {
  return out << "{o:" << line.origin << ", v:" << line.vel << "}";
}

Line Line::operator+(Vec4 delta) { return Line(origin + delta, vel); }

Line Line::operator-(Vec4 delta) { return Line(origin - delta, vel); }

Vec4 Line::PosAfter(float delta_time) const {
  return Vec4(origin.t + delta_time, origin.r + vel * delta_time);
}

Vec4 Line::PosAt(float time) const { return PosAfter(time - origin.t); }

Line Line::TransformedToFrame(const Vec3 &other_vel) const {
  return Line(origin.TransformedToFrame(other_vel),
              vel.VelTransformedToFrame(other_vel));
}

Line Line::TransformedFromFrame(const Vec3 &other_vel) const {
  return TransformedToFrame(-other_vel);
}

Line Line::TransformedToFrame(const Line &other) const {
  return Line(origin.TransformedToFrame(other),
              vel.VelTransformedToFrame(other.vel));
}

Line Line::TransformedFromFrame(const Line &other) const {
  return Line(origin.TransformedFromFrame(other),
              vel.VelTransformedFromFrame(other.vel));
}

/*
  RANDOM MATH
*/

std::mt19937 r_gen;
std::uniform_real_distribution<float> rand_real(0.0f, 1.0f);

float RandomReal() { return rand_real(r_gen); }

Vec3 RandomVectorInUnitBall() {
  Vec3 vec(1, 1, 1);
  std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
  while (vec.NormSq() >= 1 || vec.NormSq() == 0.0f) {
    vec =
        Vec3(2 * RandomReal() - 1, 2 * RandomReal() - 1, 2 * RandomReal() - 1);
  }
  return vec;
}

Vec3 RandomUnitVector() {
  Vec3 vec = RandomVectorInUnitBall();
  return vec.NormalizedNonzero();
}

Vec2 RandomVectorInUnitDisk() {
  float x = 1.0f;
  float y = 1.0f;
  while (x * x + y * y >= 1.0f) {
    x = RandomReal() * 2.0f - 1.0f;
    y = RandomReal() * 2.0f - 1.0f;
  }
  return Vec2(x, y);
}
