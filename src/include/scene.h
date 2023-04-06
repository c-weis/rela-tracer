// Copyright 2023 Christoph Weis
#ifndef SCENE_H_INCLUDED
#define SCENE_H_INCLUDED

#include <memory>
#include <optional>
#include <random>
#include <vector>

#include "materials.h"
#include "math.h"

class Object;
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

typedef std::shared_ptr<Material> MaterialPtr;
class Object {
 public:
  Object(Material mat, Line worldline)
      : worldline_(worldline), mat_(std::make_shared<Material>(mat)) {}
  Object(MaterialPtr mat, Line worldline)
      : worldline_(worldline), mat_(mat) {}

  Vec4 PosAt(float time) const;
  Vec4 PosAfter(float d_time) const;
  OptionalHitRecord intersect(const Line &) const;

  MaterialPtr GetMaterial() const { return mat_; }
  void SetMaterial(Material mat) { mat_ = std::make_shared<Material>(mat); }

  Vec4 GetOrigin() const { return worldline_.origin; }
  Vec3 GetVelocity() const { return worldline_.vel; }
  Line GetWorldline() const { return worldline_; }

 private:
  Line worldline_;
  MaterialPtr mat_;

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
typedef std::shared_ptr<Object> ObjectPointer;
// Polymorphic container for objects
typedef std::vector<ObjectPointer> ObjectList;

class Sphere : public Object {
 public:
  // Construct sphere from centre worldline and radius
  Sphere(Material mat, Line worldline, float rad)
      : Object(mat, worldline), rad_(rad) {}
  // Construct sphere from centre point `O`, radius `rad`,
  // timeslice `t`, and velocity `vel` (default: stationary)
  Sphere(Material mat, Vec3 O, float rad, float t = 0, Vec3 vel = kZero3)
      : Object(mat, Line(Vec4(t, O), vel)), rad_(rad) {}
  float rad();

 private:
  const float rad_;  // radius
  OptionalReferenceFrameHit intersect_in_rest_frame(const Line &) const;
};

class Shape2D : public Object {
 public:
  Shape2D(Material mat, Line worldline, Vec3 a, Vec3 b)
      : Object(mat, worldline), a_(a), b_(b), normal_(NormalTo(a, b)) {}

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
  Triangle(Material mat, Line worldline, Vec3 a, Vec3 b)
      : Shape2D(mat, worldline, a, b) {}
  // Construct triangle given the three vertices `A`, `B`, `C`,
  // timeslice `t` and velocity `vel` (default: stationary)
  Triangle(Material mat, Vec3 A, Vec3 B, Vec3 C, float t = 0, Vec3 vel = kZero3)
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
  Parallelogram(Material mat, Line worldline, Vec3 a, Vec3 b)
      : Shape2D(mat, worldline, a, b) {}
  // Construct Parallelogram from three adjacent vertices 'A', 'B', 'C',
  // timeslice `t` and velocity `vel` (default: stationary)
  Parallelogram(Material mat, Vec3 A, Vec3 B, Vec3 C, float t = 0,
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
  Plane(Material mat, Line worldline, Vec3 a, Vec3 b)
      : Shape2D(mat, worldline, a, b) {}
  // Construct plane by giving 3 points `A`, `B`, `C` on it,
  // a timeslice `t` and a velocity `vel`
  // Defaults to a stationary plane if `t` and `vel` are
  // not specified.
  Plane(Material mat, Vec3 A, Vec3 B, Vec3 C, float t = 0, Vec3 vel = kZero3)
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
class Box {
 private:
  Line worldline_BUR_;

 public:
  // we say `up` & `down` to avoid the Bottom-Back intials clash
  Parallelogram front;
  Parallelogram down;
  Parallelogram left;
  Parallelogram back;
  Parallelogram right;
  Parallelogram up;

  // Construct Box from worldline of Front Down Left corner
  // and edge vectors front_to_back, down_to_up, left_to_right
  Box(Material mat, Line worldline_FDL, Vec3 front_to_back, Vec3 down_to_up,
      Vec3 left_to_right)
      : worldline_BUR_(Vec4(worldline_FDL.origin.t,
                            worldline_FDL.origin.r + front_to_back +
                                down_to_up + left_to_right),
                       worldline_FDL.vel),
        front(mat, worldline_FDL, left_to_right, down_to_up),
        down(mat, worldline_FDL, front_to_back, left_to_right),
        left(mat, worldline_FDL, down_to_up, front_to_back),
        back(mat, worldline_BUR_, -down_to_up, -left_to_right),
        up(mat, worldline_BUR_, -left_to_right, -front_to_back),
        right(mat, worldline_BUR_, -front_to_back, -down_to_up) {}

  // Construct box given 4 vertices `FDL` (front, down, left),
  // `BDL` (back, down, left), `FUL` (front, up, left),
  // `FDR` (front, down, right) and optionally
  // a timeslice `t` and velocity `vel`. (Default=stationary)
  Box(Material mat, Vec3 FDL, Vec3 BDL, Vec3 FUL, Vec3 FDR, float t = 0,
      Vec3 vel = kZero3)
      : Box(mat, Line(Vec4(t, FDL), vel), BDL - FDL, FUL - FDL, FDR - FDL) {}
};

class Camera {
 public:
  // Constructs camera from `worldline` of camera sensor,
  // and screen data (given in the restframe of the camera
  // with the camera at the origin).
  Camera(Line worldline, Vec3 cam_to_screen_centre, Vec3 screen_centre_to_right,
         int width, int height)
      : worldline_(worldline),
        screen_centre_(cam_to_screen_centre),
        screen_right_(screen_centre_to_right),
        screen_width_(width),
        screen_height_(height) {}

  Camera(Vec3 origin, Vec3 screen_centre, Vec3 screen_centre_to_right,
         int width, int height, float t = 0, Vec3 vel = kZero3)
      : Camera(Line(Vec4(t, origin), vel), screen_centre,
               screen_centre_to_right, width, height) {}

  int GetScreenWidth() const { return screen_width_; }
  int GetScreenHeight() const { return screen_height_; }
  Vec3 GetVelocity() const { return worldline_.vel; }
  // Get image rays at specified point `t_cam_frame` in time
  // (given in the restframe of the camera, defaults to 0)
  LineList ImageRays(int rays_per_pixel, std::mt19937 *r_gen,
                     float t_cam_frame = 0) const;

 private:
  // worldline of (pointlike) camera sensor
  Line worldline_;

  // screen data
  Vec3 screen_centre_;
  Vec3 screen_right_;  // middle of right screen edge
  // screen size in pixels
  int screen_width_;
  int screen_height_;
};
typedef std::vector<Camera> CameraList;

class Scene {
 public:
  template <class Derived_Object>
  void Add(const Derived_Object &_obj) {
    objects_.push_back(std::make_shared<Derived_Object>(_obj));
  }

  void Add(const Box &box) {
    Add(box.front);
    Add(box.down);
    Add(box.left);
    Add(box.back);
    Add(box.up);
    Add(box.right);
  }

  void AddCamera(const Camera &);

  void SetSkylight(Material material, Vec3 skylight_direction) {
    skylight_ = std::make_shared<Material>(material);
    skylight_direction_ = skylight_direction;
  }

  const Camera &GetCamera(int index = 0) const;

  // Finds the nearest (most recent) intersection of ray with any object
  // in the scene
  OptionalHitRecord MostRecentHit(const Line &ray) const;

  ColorData BackgroundColor(const Line &ray) const;

 private:
  CameraList cameras_;
  ObjectList objects_;
  MaterialPtr skylight_ = nullptr;
  Vec3 skylight_direction_ = kDefaultNormal;
};

#endif