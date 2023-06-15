// Copyright 2023 Christoph Weis
#ifndef MATH_H_INCLUDED
#define MATH_H_INCLUDED

/*
  BASIC 3- AND 4-VECTOR MATH
*/

#include <iostream>
#include <optional>
#include <vector>

// Forward declarations
struct Vec3;
struct Vec4;
struct Line;
typedef std::optional<Vec3> OptionalVec3;
typedef std::optional<Vec4> OptionalVec4;

// Avoid division by tiny numbers by using below.
const float kDivisionEpsilon = 1e-20;

// 3-vector class implemented with a view towards special relativity
struct Vec3 {
  float x, y, z;
  Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
  friend std::ostream &operator<<(std::ostream &, const Vec3 &);

  Vec3 NormalizedNonzero() const;
  OptionalVec3 Normalized() const;

  float NormSq() const;
  float Norm() const;
  float Dot3(const Vec3 &) const;
  float CosAngleTo(const Vec3 &) const;
  Vec3 Cross(const Vec3 &) const;

  void NormalizeNonzero();
  bool Normalize();

  // calculate Gamma factor of 3-vector considered as velocity
  float Gamma() const;
  // calculate 1/Gamma factor of 3-vector considered as velocity
  float GammaInv() const;

  // transform 3-vector (considered as a velocity) to frame with given relative
  // velocity
  Vec3 VelTransformedToFrame(const Vec3 &rel_vel) const;
  // transform 3-vector (considered as a velocity) from frame with given
  // relative velocity
  Vec3 VelTransformedFromFrame(const Vec3 &rel_vel) const;

  // reflect 3-vector across plane with given normal
  Vec3 Reflect(const Vec3 &plane_normal) const;
};
Vec3 operator+(const Vec3 &, const Vec3 &);
Vec3 operator-(const Vec3 &, const Vec3 &);
Vec3 operator-(const Vec3 &);
Vec3 operator*(const Vec3 &, float);
Vec3 operator*(float, const Vec3 &);
Vec3 operator/(const Vec3 &, float);
float Dot3(const Vec3 &, const Vec3 &);
// Cross product
Vec3 operator^(const Vec3 &, const Vec3 &);
Vec3 Cross(const Vec3 &left, const Vec3 &right);
Vec3 NormalTo(const Vec3 &left, const Vec3 &right);

const Vec3 kZero3(0, 0, 0);
// Hack: default normal vector when trying to compute
//       normal to colinear pair of vectors
const Vec3 kDefaultNormal = Vec3(1, 0, 0);

// 4-vector class, implemented with view towards special relativity
struct Vec4 {
  float t;  // time coordinate
  Vec3 r;   // space coordinate
  Vec4(float _t, float _x, float _y, float _z) : t(_t), r(Vec3(_x, _y, _z)) {}
  Vec4(float _t, const Vec3 &_s) : t(_t), r(_s) {}
  friend std::ostream &operator<<(std::ostream &, const Vec4 &);

  float Dot4(const Vec4 &) const;

  // Returns -t^2 + r^2.
  float MinkowskiNormSq() const;
  // Returns r^2.
  float SpaceNormSq() const;
  // Returns sqrt(-t^2 + r^2).
  float MinkowskiNorm() const;
  // Returns |r|
  float SpaceNorm() const;

  // Boost 4-vector to frame with given relative velocity.
  Vec4 TransformedToFrame(const Vec3 &rel_vel) const;  // boost
  // De-boost 4-vector from frame with given relative velocity.
  Vec4 TransformedFromFrame(const Vec3 &rel_vel) const;  // de-boost
  // Transform 4-vector to frame with given worldline (subtracts origin of
  // `rf_line` before transforming relative 4-vector).
  Vec4 TransformedToFrame(const Line &rf_line) const;  // boost
  // Transform 4-vector from frame with given worldline (adds origin of
  // `rf_line` after transforming 4-vector).
  Vec4 TransformedFromFrame(const Line &rf_line) const;  // de-boost
};
Vec4 operator+(const Vec4 &, const Vec4 &);
Vec4 operator-(const Vec4 &, const Vec4 &);
Vec4 operator*(const Vec4 &, float);
Vec4 operator/(const Vec4 &, float);
Vec4 operator-(const Vec4 &);
Vec4 operator*(float, const Vec4 &);
float Dot4(const Vec4 &, const Vec4 &);

const Vec4 kZero4(0, kZero3);

// Line through spacetime.
struct Line {
  Vec4 origin;  // ray origin (in spacetime)
  Vec3 vel;     // velocity

  Line(Vec4 _origin, Vec3 _vel) : origin(_origin), vel(_vel) {}
  friend std::ostream &operator<<(std::ostream &, const Line &);

  Line operator+(Vec4 delta);
  Line operator-(Vec4 delta);

  // Returns origin + delta_time * (1, vel)
  Vec4 PosAfter(float delta_time) const;
  // Returns 4-position of this line with given time coordinate.
  Vec4 PosAt(float time) const;

  // Boosts line to frame with given relative velocity.
  Line TransformedToFrame(const Vec3 &rel_vel) const;  // boost
  // De-boosts line from frame with given relative velocity.
  Line TransformedFromFrame(const Vec3 &rel_vel) const;  // de-boost
  // Transforms line to frame with given worldline (subtracts origin
  // before transforming relative 4-vector).
  Line TransformedToFrame(const Line &rf_line) const;  // boost
  // Transforms line to frame with given worldline (adds origin
  // after transforming 4-vector).
  Line TransformedFromFrame(const Line &rf_line) const;  // de-boost
};
typedef std::vector<Line> LineList;

// Struct containing restframe data about a hit.
struct ReferenceFrameHit {
  Vec4 pos;        // 4-position of hit
  Vec3 normal;     // normal of object at hit point
  Vec3 scattered;  // velocity of scattered ray
  bool outside;    // bool specifying whether we hit the outside of the object

  ReferenceFrameHit(Vec4 _pos, Vec3 _normal, Vec3 _scattered,
                    bool _outside = true)
      : pos(_pos), normal(_normal), scattered(_scattered), outside(_outside) {}
};
typedef std::optional<ReferenceFrameHit> OptionalReferenceFrameHit;

// Forward declaration needed for HitRecord.
class Object;
// Record of a ray hitting an object. Contains a pointer to the hit object and
// hit data in the rest frame of the object.
struct HitRecord {
  // Constructor
  // Args:
  //  rf_pos: position of hit in 4-space, given in object restframe
  //  rf_normal: object normal at hit position, given in object restframe
  //  rf_scattered: the velocity of the ray emitted by the ray, given in object
  //   restframe
  //  obj: pointer to object hit
  HitRecord(Vec4 rf_pos, Vec3 rf_normal, Vec3 rf_scattered, const Object *obj)
      : obj(obj), rf(rf_pos, rf_normal, rf_scattered) {}
  // Constructor
  // Args:
  //  rf_pos: position of hit in 4-space, given in object restframe
  //  rf_normal: object normal at hit position, given in object restframe
  //  rf_scattered: the velocity of the ray emitted by the ray, given in object
  //   restframe
  //  obj: pointer to object hit
  HitRecord(ReferenceFrameHit rf, const Object *obj) : obj(obj), rf(rf) {}
  // Empty destructor: *obj is not owned by HitRecord
  ~HitRecord() = default;
  // Default copy constructor.
  HitRecord(const HitRecord &other) = default;
  // Default assignment operator.
  HitRecord &operator=(const HitRecord &) = default;

  // time passed / ray length in standard(!) frame
  float HitTime() const;

  // Pointer to the hit object.
  const Object *obj;
  // HitRecord record in restframe of *obj_.
  ReferenceFrameHit rf;
};
typedef std::optional<HitRecord> OptionalHitRecord;

// RANDOMNESS MATH

// Uniformly draws float float in [0,1).
float RandomReal();
// Uniformly draws vector from unit 2-sphere.
Vec3 RandomUnitVector();
// Uniformly draws vector from unit 3-ball.
Vec3 RandomVectorInUnitBall();

struct Vec2 {
  Vec2(float _x, float _y) : x(_x), y(_y) {}
  float x, y;
};

// Uniformly draws vector from unit 2-disk.
Vec2 RandomVectorInUnitDisk();

#endif
