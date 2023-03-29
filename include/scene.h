// Copyright 2023 Christoph Weis
#pragma once

#include <memory>
#include <optional>
#include <random>
#include <vector>

#include "include/math.h"

class Object;
class Material;

struct HitRecord {
  // Pointer to the hit object, *obj is not owned by HitRecord.
  const Object *obj;
  // HitRecord record in restframe of *obj_.
  ReferenceFrameHit rf;

  HitRecord(Vec4 rf_pos, Vec3 rf_normal, Vec3 rf_scattered, const Object *obj)
      : obj(obj), rf(rf_pos, rf_normal, rf_scattered) {}
  HitRecord(ReferenceFrameHit rf, const Object *obj) : obj(obj), rf(rf) {}

  ~HitRecord() =
      default;  // empty destructor because *obj is not owned by HitRecord
  HitRecord(const HitRecord &other) = default;
  HitRecord &operator=(const HitRecord &) = default;

  // time passed / ray length in !standard! frame
  float hit_time() const;
};
typedef std::optional<HitRecord> OptionalHitRecord;

class Object {
 public:
  const Line worldline_;
  const Material *mat_;

  Object(Material *mat, Line worldline) : worldline_(worldline), mat_(mat) {}

  OptionalVec4 PosAt(float time) const;
  OptionalVec4 PosAfter(float d_time) const;
  OptionalHitRecord intersect(const Line &) const;

 private:
  virtual OptionalReferenceFrameHit intersect_in_rest_frame(
      const Line &) const = 0;
};

struct UVCoordinates {
  float u;
  float v;
};

// Return t such that x_0 + vel * t intersects
// the plane through the origin with given normal
std::optional<float> plane_intersection_time(Vec3 x_0, Vec3 vel,
                                             Vec3 plane_normal);

// Given vectors x, a, b, returns u,v such that
// x' = u * a + v * b is the closest approximation to x
// in the plane spanned by a and b
// (mostly used when x is on the plane spanned by a and b)
UVCoordinates uv_coordinates(Vec3 x, Vec3 a, Vec3 b);

// Return ReferenceFrameHit for think object, automating
// INSIDE/OUTSIDE detection
ReferenceFrameHit rf_hit_from(Vec4 pos, Vec3 ray_vel, Vec3 outside_normal);

// Polymorphic pointer
typedef std::unique_ptr<Object> ObjectPointer;
// Polymorphic container for objects
typedef std::vector<std::unique_ptr<Object>> ObjectList;

class Sphere : public Object {
 public:
  // Construct sphere from centre worldline and radius
  Sphere(Material *mat, Line worldline, float rad)
      : Object(mat, worldline), rad_(rad) {}
  // Construct sphere from centre point `O`, radius `rad`,
  // timeslice `t`, and velocity `vel` (default: stationary)
  Sphere(Material *mat, Vec3 O, float rad, float t = 0, Vec3 vel = kZero3)
      : Object(mat, Line(Vec4(t, O), vel)), rad_(rad) {}
  float rad();

 private:
  const float rad_;  // radius
  OptionalReferenceFrameHit intersect_in_rest_frame(const Line &) const;
};

class Shape2D : public Object {
 public:
  Shape2D(Material *mat, Line worldline, Vec3 a, Vec3 b)
      : a_(a), b_(b), normal_(NormalTo(a, b)), Object(mat, worldline) {}

 private:
  const Vec3 a_, b_;   // two vectors spanning the shape
  const Vec3 normal_;  // outside normal of the shape

  OptionalReferenceFrameHit intersect_in_rest_frame(const Line &) const;
  bool is_inside_shape(UVCoordinates uv) const {
    return is_inside_shape(uv.u, uv.v);
  }

  virtual bool is_inside_shape(float u, float v) const = 0;
};

class Triangle : public Shape2D {
 public:
  // Construct triangle from `worldline` of one vertex
  // and sides `a`, `b` in restframe of that vertex
  Triangle(Material *mat, Line worldline, Vec3 a, Vec3 b)
      : Shape2D(mat, worldline, a, b) {}
  // Construct triangle given the three vertices `A`, `B`, `C`,
  // timeslice `t` and velocity `vel` (default: stationary)
  Triangle(Material *mat, Vec3 A, Vec3 B, Vec3 C, float t = 0,
           Vec3 vel = kZero3)
      : Shape2D(mat, Line(Vec4(t, A), vel), B - A, C - A) {}
  // The worldline of the shape tracks one of the vertices.
  // a_ and b_ are the positions of the other vertices in
  // the restframe of the triangle.
 private:
  // check if u > 0, v > 0, u+v <= 1
  bool is_inside_shape(float u, float v) const override;
};

class Parallelogram : public Shape2D {
 public:
  // Construct Parallelogram from worldline of a vertex 0
  // and `a`, `b`, the two adjacent vertices of 0
  // in the restframe (with 0 at the origin)
  Parallelogram(Material *mat, Line worldline, Vec3 a, Vec3 b)
      : Shape2D(mat, worldline, a, b) {}
  // Construct Parallelogram from three adjacent vertices 'A', 'B', 'C',
  // timeslice `t` and velocity `vel` (default: stationary)
  Parallelogram(Material *mat, Vec3 A, Vec3 B, Vec3 C, float t = 0,
                Vec3 vel = kZero3)
      : Shape2D(mat, Line(Vec4(t, A), vel), B - A, C - A) {}

 private:
  // check if u > 0, v > 0, u < 1, v < 1
  bool is_inside_shape(float u, float v) const override;
};

class Plane : public Shape2D {
 public:
  // Construct plane from worldline of a point
  // and two vectors `a`, `b` spanning the plane
  // in the restframe of that point
  Plane(Material *mat, Line worldline, Vec3 a, Vec3 b)
      : Shape2D(mat, worldline, a, b) {}
  // Construct plane by giving 3 points `A`, `B`, `C` on it,
  // a timeslice `t` and a velocity `vel`
  // Defaults to a stationary plane if `t` and `vel` are
  // not specified.
  Plane(Material *mat, Vec3 A, Vec3 B, Vec3 C, float t = 0, Vec3 vel = kZero3)
      : Shape2D(mat, Line(Vec4(t, A), vel), B - A, C - A){};

 private:
  bool is_inside_shape(float u, float v) const override;
};

// TODO(c): implement composites
// An object comprised of multiple sub-objects,
// all moving uniformly - simplify intersection logic by
// turning Lorentz transformations into simple shifts
// class CompositeObject : public Object {
//};

class Camera {
 public:
  Camera(Vec3 _cam_pos, float _cam_time, Vec3 _screen_centre,
         Vec3 _screen_right, int _width, int _height)
      : cam_pos_(_cam_pos),
        cam_time_(_cam_time),
        screen_centre_(_screen_centre),
        screen_right_(_screen_right),
        screen_width_(_width),
        screen_height_(_height) {}

  Camera(Vec4 _cam_pos4, Vec3 _screen_centre, Vec3 _screen_right, int _width,
         int _height)
      : cam_pos_(_cam_pos4.r),
        cam_time_(_cam_pos4.t),
        screen_centre_(_screen_centre),
        screen_right_(_screen_right),
        screen_width_(_width),
        screen_height_(_height) {}

  int ScreenWidth() const;
  int ScreenHeight() const;
  LineList ImageRays(int rays_per_pixel, std::mt19937 *r_gen) const;

 private:
  // camera 4-vector, separated for convenience
  Vec3 cam_pos_;
  float cam_time_;

  // screen data
  Vec3 screen_centre_;
  Vec3 screen_right_;  // middle of right screen edge
  // screen size in pixels
  int screen_width_;
  int screen_height_;
};
typedef std::vector<Camera> CameraList;

struct Scene {
 public:
  template <class Derived_Object>
  void Add(const Derived_Object &_obj) {
    objects_.push_back(std::make_unique<Derived_Object>(_obj));
  }
  void AddCamera(const Camera &);

  // Adds a box to the scene (as a collection of parallelograms)
  // given the `worldline` of one vertex and the sides `a`, `b`, `c`
  // of the box
  void AddBox(Material *mat, Line worldline, Vec3 a, Vec3 b, Vec3 c);
  // Add box by specifying 4 vertices `O` (origin), `A`, `B`, `C`,
  // timeslice `t` and velocity `vel (default: stationary)
  void AddBox(Material *mat, Vec3 O, Vec3 A, Vec3 B, Vec3 C, float t = 0,
              Vec3 vel = kZero3);
  const Camera &GetCamera(int index = 0) const;

  OptionalHitRecord MostRecentHit(const Line &ray) const;
  /*
  LineList scattered(const Line &ray)
      const;  // takes a Line, finds nearest intersection, scatters
  */

 private:
  CameraList cameras_;
  ObjectList objects_;
};
