// Copyright 2023 Christoph Weis
#ifndef MATH_H_INCLUDED
#define MATH_H_INCLUDED

#include <iostream>
#include <optional>
#include <vector>

struct Vec3;
struct Vec4;
struct Line;
typedef std::optional<Vec3> OptionalVec3;
typedef std::optional<Vec4> OptionalVec4;

// Avoid division by tiny numbers by using below.
const float kDivisionEpsilon = 1e-20;

struct Vec3 {
  float x, y, z;
  Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
  friend std::ostream& operator<< (std::ostream&, const Vec3 &);

  Vec3 NormalizedNonzero() const;
  OptionalVec3 Normalized() const;

  float NormSq() const;
  float Norm() const;
  float Dot3(const Vec3 &) const;
  float CosAngleTo(const Vec3 &) const;
  Vec3 Cross(const Vec3 &) const;

  void NormalizeNonzero();
  bool Normalize();

  // considered as relative velocities
  float GammaInv() const;
  float Gamma() const;

  Vec3 TransformedToFrame(const Vec3 &) const;
  Vec3 TransformedFromFrame(const Vec3 &) const;
};
Vec3 operator+(const Vec3 &, const Vec3 &);
Vec3 operator-(const Vec3 &, const Vec3 &);
Vec3 operator-(const Vec3 &);
Vec3 operator*(const Vec3 &, float);
Vec3 operator*(float, const Vec3 &);
Vec3 operator/(const Vec3 &, float);
float Dot3(const Vec3 &, const Vec3 &);  // dot product
Vec3 operator^(const Vec3 &, const Vec3 &);  // cross product 
Vec3 Cross(const Vec3 &left, const Vec3 &right);
Vec3 NormalTo(const Vec3 &left, const Vec3 &right);

const Vec3 kZero3(0,0,0);
// Hack: default normal vector when tryuing to compute
//       normal to colinear pair of vectors
const Vec3 kDefaultNormal = Vec3(1, 0, 0);  // bad hack

struct Vec4 {
  float t;  // time coordinate
  Vec3 r;   // space coordinate
  Vec4(float _t, float _x, float _y, float _z) : t(_t), r(Vec3(_x, _y, _z)) {}
  Vec4(float _t, const Vec3 &_s) : t(_t), r(_s) {}
  friend std::ostream& operator<< (std::ostream&, const Vec4 &);

  float Dot4(const Vec4 &) const;

  float MinkowskiNormSq() const;
  float SpaceNormSq() const;
  float MinkowskiNorm() const;
  float SpaceNorm() const;

  Vec4 TransformedToFrame(const Vec3 &) const;    // boost
  Vec4 TransformedFromFrame(const Vec3 &) const;  // de-boost
  Vec4 TransformedToFrame(const Line &) const;    // boost
  Vec4 TransformedFromFrame(const Line &) const;  // de-boost
};
Vec4 operator+(const Vec4 &, const Vec4 &);
Vec4 operator-(const Vec4 &, const Vec4 &);
Vec4 operator*(const Vec4 &, float);
Vec4 operator/(const Vec4 &, float);
Vec4 operator-(const Vec4 &);
Vec4 operator*(float, const Vec4 &);
float Dot4(const Vec4 &, const Vec4 &);

const Vec4 kZero4(0, kZero3);

/*
    Line through spacetime
*/

struct Line {
  Vec4 origin;    // ray origin (in spacetime)
  Vec3 vel;       // velocity

  Line(Vec4 _origin, Vec3 _vel)
      : origin(_origin), vel(_vel) {}
  friend std::ostream& operator<< (std::ostream&, const Line &);

  Vec4 PosAfter(float delta_time) const;
  Vec4 PosAt(float time) const;

  // Lorenz boosts
  Line TransformedToFrame(const Vec3 &) const;    // boost
  Line TransformedFromFrame(const Vec3 &) const;  // de-boost
  Line TransformedToFrame(const Line &) const;    // boost
  Line TransformedFromFrame(const Line &) const;  // de-boost
};
typedef std::vector<Line> LineList;

struct ReferenceFrameHit {
  /*
      A struct for restframe hit data.
  */
  Vec4 pos;
  Vec3 normal;
  Vec3 scattered;  // velocity of scattered ray
  bool outside;  // did we hit the outside or inside?

  ReferenceFrameHit(Vec4 _pos, Vec3 _normal, Vec3 _scattered, 
  bool _outside=true)
      : pos(_pos), normal(_normal), scattered(_scattered), outside(_outside) {}
};
typedef std::optional<ReferenceFrameHit> OptionalReferenceFrameHit;

#endif