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
  float HitTime() const;
};
typedef std::optional<HitRecord> OptionalHitRecord;

class Object {
 public:
  Object(Line worldline, Material mat = Material())
      : worldline_(worldline), mat_(mat) {}

  Vec4 PosAt(float time) const;
  Vec4 PosAfter(float d_time) const;
  OptionalHitRecord intersect(const Line &) const;

  const Material *GetMaterial() const { return &mat_; }
  void SetMaterial(Material mat) { mat_ = mat; }

  Vec4 GetOrigin() const { return worldline_.origin; }
  Vec3 GetVelocity() const { return worldline_.vel; }
  Line GetWorldline() const { return worldline_; }

 private:
  Line worldline_;
  Material mat_;

  virtual OptionalReferenceFrameHit intersect_in_rest_frame(
      const Line &) const = 0;
};

struct UVCoordinates {
  float u;
  float v;
};

// Return t such that x_0 + vel * t intersects
// the plane through the origin with given normal
std::optional<float> plane_intersection_time(Vec3 x_0, Vec3 vel, Vec3 plane_o,
                                             Vec3 plane_normal);

// Given vectors x, o, a, b, returns u,v such that
// x' = u * a + v * b is the closest approximation to x
// in the plane going through o spanned by a and b
// (mostly used when x is on the plane spanned by a and b)
UVCoordinates uv_coordinates(Vec3 x, Vec3 o, Vec3 a, Vec3 b);

// Return ReferenceFrameHit for thin object, automating
// INSIDE/OUTSIDE detection
ReferenceFrameHit rf_hit_from(Vec4 pos, Vec3 ray_vel, Vec3 outside_normal);

// Polymorphic pointer
typedef std::shared_ptr<Object> ObjectPointer;
// Polymorphic container for objects
typedef std::vector<ObjectPointer> ObjectList;

class Sphere : public Object {
 public:
  // Construct sphere from center worldline and restframe radius
  Sphere(Line worldline, float rad, Material mat = Material())
      : Object(worldline, mat), rad_(rad) {}
  // Construct sphere from center point `O`, restframe radius `rad`,
  // timeslice `t`, and velocity `vel` (default: stationary)
  Sphere(Vec3 O, float rad, Vec3 vel = kZero3, float t = 0,
         Material mat = Material())
      : Object(Line(Vec4(t, O), vel), mat), rad_(rad) {}
  Sphere(Vec3 O, float rad, Material mat) : Sphere(O, rad, kZero3, 0, mat) {}
  float rad();

 private:
  const float rad_;  // radius
  OptionalReferenceFrameHit intersect_in_rest_frame(const Line &) const;
};

class Shape2D : public Object {
 public:
  Shape2D(Line worldline, Vec3 o, Vec3 a, Vec3 b, Material mat = Material())
      : Object(worldline, mat), o_(o), a_(a), b_(b), normal_(NormalTo(a, b)) {}

 private:
  // all in rest frame:
  const Vec3 o_;       // origin of shape in rest frame
  const Vec3 a_, b_;   // two vectors spanning the plane of the shape
  const Vec3 normal_;  // outside normal of the shape, (a ∧ b).Normalized()

  OptionalReferenceFrameHit intersect_in_rest_frame(const Line &) const;
  bool is_inside_shape(UVCoordinates uv) const {
    return is_inside_shape(uv.u, uv.v);
  }

  virtual bool is_inside_shape(float u, float v) const = 0;
};

class Triangle : public Shape2D {
 public:
  // Construct triangle from `worldline` of a point
  // and the three vertices `A_rf`, `B_rf`, `C_rf` in restframe
  Triangle(Line worldline, Vec3 A_rf, Vec3 B_rf, Vec3 C_rf,
           Material mat = Material())
      : Shape2D(worldline, A_rf, B_rf - A_rf, C_rf - A_rf, mat) {}
  // Construct triangle given a point `O` on it in the standard frame,
  // the three vertices `A_rf`, `B_rf`, `C_rf` in its restframe (relative to O),
  // velocity `vel` of O (default: stationary)
  // and the time `t` at which its position is given (default: 0)
  Triangle(Vec3 O, Vec3 A_rf, Vec3 B_rf, Vec3 C_rf, Vec3 vel = kZero3,
           float t = 0, Material mat = Material())
      : Triangle(Line(Vec4(t, O), vel), A_rf, B_rf, C_rf, mat) {}

 private:
  // check if u > 0, v > 0, u+v <= 1
  bool is_inside_shape(float u, float v) const override;
};

class Parallelogram : public Shape2D {
 public:
  // Construct Parallelogram from the `worldline` of the coordinate origin,
  // the position `o_rf` of a vertex in its restframe,
  // and `a_rf`, `b_rf`, its spanning sides in the restframe
  Parallelogram(Line worldline, Vec3 o_rf, Vec3 a_rf, Vec3 b_rf,
                Material mat = Material())
      : Shape2D(worldline, o_rf, a_rf, b_rf, mat) {}
  // Construct Parallelogram from `worldline_c` of the center point,
  // vectors `half_a`, `half_b` from center to the midpoints of
  // the sides of the parallelogram
  Parallelogram(Line worldline_c, Vec3 half_a, Vec3 half_b,
                Material mat = Material())
      : Parallelogram(worldline_c, -(half_a + half_b), 2 * half_a, 2 * half_b,
                      mat) {}
  // Construct Parallelogram from position `C` of the center point,
  // vectors `half_a`, `half_b` from center to the midpoints of
  // the sides of the parallelogram, velocity `vel` (default: stationary)
  // and timeslice `t` (default: 0)
  Parallelogram(Vec3 C, Vec3 half_a, Vec3 half_b, Vec3 vel = kZero3,
                float t = 0.0f, Material mat = Material())
      : Parallelogram(Line(Vec4(t, C), vel), -(half_a + half_b), 2 * half_a,
                      2 * half_b, mat) {}
  Parallelogram(Vec3 C, Vec3 half_a, Vec3 half_b, Material mat = Material())
      : Parallelogram(Line(Vec4(0, C), kZero3), -(half_a + half_b), 2 * half_a,
                      2 * half_b, mat) {}

 private:
  // check if u > 0, v > 0, u < 1, v < 1
  bool is_inside_shape(float u, float v) const override;
};

class Plane : public Shape2D {
 public:
  // Construct plane from worldline of a point
  // and two vectors `a_rf`, `b_rf` spanning the plane
  // in the restframe of that point
  Plane(Line worldline, Vec3 a_rf, Vec3 b_rf, Material mat = Material())
      : Shape2D(worldline, kZero3, a_rf, b_rf, mat) {}
  // Construct plane by giving a point `O`, on it (in the standard frame),
  // two vectors `a_rf`, `b_rf` spanning it in its restframe,
  // the velocity `vel` (default: stationary) and timeslice `t` (default: 0)
  Plane(Vec3 O, Vec3 a_rf, Vec3 b_rf, Vec3 vel = kZero3, float t = 0,
        Material mat = Material())
      : Plane(Line(Vec4(t, O), vel), a_rf, b_rf, mat) {}
  Plane(Vec3 O, Vec3 a_rf, Vec3 b_rf, Material mat = Material())
      : Plane(O, a_rf, b_rf, kZero3, 0, mat) {}

 private:
  bool is_inside_shape(float u, float v) const override;
};

class Box {
 public:
  // we say `up` & `down` to avoid the Bottom-Back intials clash
  Parallelogram right;
  Parallelogram back;
  Parallelogram up;
  Parallelogram left;
  Parallelogram front;
  Parallelogram down;

  // Construct Box from `worldline_c` of center
  // and vectors `c_to_right`, `c_to_back`, `c_to_up`
  // from the center to the midpoints of right, back
  // and up face, respectively.
  Box(Line worldline_c, Vec3 c_to_right, Vec3 c_to_back, Vec3 c_to_up,
      Material mat = Material())
      : right(worldline_c, c_to_right + c_to_back + c_to_up, -2 * c_to_back,
              -2 * c_to_up, mat),
        back(worldline_c, c_to_right + c_to_back + c_to_up, -2 * c_to_up,
             -2 * c_to_right, mat),
        up(worldline_c, c_to_right + c_to_back + c_to_up, -2 * c_to_right,
           -2 * c_to_back, mat),
        left(worldline_c, -c_to_right - c_to_back - c_to_up, 2 * c_to_up,
             2 * c_to_back, mat),
        front(worldline_c, -c_to_right - c_to_back - c_to_up, 2 * c_to_right,
              2 * c_to_up, mat),
        down(worldline_c, -c_to_right - c_to_back - c_to_up, 2 * c_to_back,
             2 * c_to_right, mat) {}

  // Construct box by giving its center point `O`, (in the standard frame),
  // vectors `c_to_right`, `c_to_back`, `c_to_up`
  // from the center to the midpoints of right, back
  // and up face, respectively.
  // Optional:
  // Velocity `vel` (default=stationary),
  // and timeslice `t` (default=0).
  Box(Vec3 O, Vec3 c_to_right, Vec3 c_to_back, Vec3 c_to_up, Vec3 vel = kZero3,
      float t = 0.0f, Material mat = Material())
      : Box(Line(Vec4(t, O), vel), c_to_right, c_to_back, c_to_up, mat) {}

  Box(Vec3 O, Vec3 c_to_right, Vec3 c_to_back, Vec3 c_to_up, Material mat)
      : Box(O, c_to_right, c_to_back, c_to_up, kZero3, 0, mat) {}
};

class Camera {
 public:
  // Constructs camera from `worldline` of camera sensor,
  // and screen data (given in the restframe of the camera
  // with the camera at the origin).
  Camera(Line worldline, Vec3 cam_to_screen_center, Vec3 screen_center_to_right,
         int width, int height)
      : worldline_(worldline),
        screen_center_(cam_to_screen_center),
        screen_center_to_right_(screen_center_to_right),
        screen_width_(width),
        screen_height_(height) {}

  Camera(Vec3 origin, Vec3 cam_to_screen_center, Vec3 screen_center_to_right,
         int width, int height, Vec3 vel = kZero3, float t = 0)
      : Camera(Line(Vec4(t, origin), vel), cam_to_screen_center,
               screen_center_to_right, width, height) {}

  int GetScreenWidth() const { return screen_width_; }
  int GetScreenHeight() const { return screen_height_; }
  Vec3 GetVelocity() const { return worldline_.vel; }
  // Get image rays at specified point `t_cam_frame` in time
  // (given in the restframe of the camera, defaults to 0)
  LineList ImageRays(float delta_t_cam_frame = 0) const;

 private:
  // worldline of (pointlike) camera sensor
  Line worldline_;

  // screen data
  Vec3 screen_center_;
  // vector from screen center to middle of right edge
  Vec3 screen_center_to_right_;
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

  Scene &AddCamera(const Camera &);

  Scene &SetSkylight(const Spectrum skylight, Vec3 skylight_direction, float skylight_exponent = 1.0f);

  const Camera &GetCamera(int index = 0) const;

  // Finds the nearest (most recent) intersection of ray with any object
  // in the scene
  OptionalHitRecord MostRecentHit(const Line &ray) const;

  Spectrum EscapedRayColor(const Line &ray) const;
  Scene &SetAmbientBackground(const Spectrum &bg_spectrum);

 private:
  CameraList cameras_;
  ObjectList objects_;
  Spectrum ambient_background_ = kBlack;
  Spectrum skylight_ = kBlack;
  Vec3 skylight_direction_ = kDefaultNormal;
  float skylight_exponent_ = 1.0f;
};

#endif