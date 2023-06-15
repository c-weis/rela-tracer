// Copyright 2023 Christoph Weis
#ifndef SCENE_H_INCLUDED
#define SCENE_H_INCLUDED

#include <memory>
#include <optional>
#include <random>
#include <vector>

#include "materials.h"
#include "math.h"

// Incomplete base-class for Objects to be placed in a scene.
// Descendants must override `intersect_in_rest_frame`.
class Object {
 public:
  // Object constructor.
  // Args:
  //  worldline: worldline of the origin of the rest frame of the object
  //  mat: Material of the object
  explicit Object(Line worldline, Material mat = Material())
      : worldline_(worldline), mat_(mat) {}

  // Returns HitRecord of closest intersection with given line
  // or std::nullopt if there is none
  OptionalHitRecord intersect(const Line &) const;

  const Material *GetMaterial() const { return &mat_; }
  void SetMaterial(Material mat) { mat_ = mat; }

  Vec4 GetOrigin() const { return worldline_.origin; }
  Vec3 GetVelocity() const { return worldline_.vel; }
  Line GetWorldline() const { return worldline_; }

 private:
  Line worldline_;
  Material mat_;

  // Intersection logic in rest frame of object. This implements the actual
  // geometry of the object. It is to be overridden by explicit deriving object
  // classes.
  virtual OptionalReferenceFrameHit intersect_in_rest_frame(
      const Line &) const = 0;
};

// Minimal UV coordinate struct.
struct UVCoordinates {
  float u;
  float v;
};

// Polymorphic pointer
typedef std::shared_ptr<Object> ObjectPointer;
// Polymorphic container for objects
typedef std::vector<ObjectPointer> ObjectList;

// Sphere specialisation of object base class.
class Sphere : public Object {
 public:
  // Sphere constructor
  // Args:
  //  worldline: worldline of center
  //  rad: sphere radius in rest frame
  //  mat: sphere material
  Sphere(Line worldline, float rad, Material mat = Material())
      : Object(worldline, mat), rad_(rad) {}
  // Sphere constructor
  // Args:
  //  O: 3-position of center (in std frame)
  //  rad: sphere radius in rest frame
  //  vel: velocity of center (in std frame)
  //  t: time at which the sphere center is at O (in std frame)
  //  mat: sphere material
  Sphere(Vec3 O, float rad, Vec3 vel = kZero3, float t = 0,
         Material mat = Material())
      : Object(Line(Vec4(t, O), vel), mat), rad_(rad) {}
  // Sphere constructor (stationary)
  // Args:
  //  O: 3-position of center (in std frame)
  //  rad: sphere radius in rest frame
  //  mat: sphere material
  Sphere(Vec3 O, float rad, Material mat) : Sphere(O, rad, kZero3, 0, mat) {}
  float rad();

 private:
  const float rad_;  // radius
  OptionalReferenceFrameHit intersect_in_rest_frame(const Line &) const;
};

// Abstract specialisation of `Object` for 2D shapes.
// Specific shapes are required to override the method `is_inside_shape`
// function, which takes two coordinates `u` and `v` and returns true if
// that point is inside the shape.
class Shape2D : public Object {
 public:
  // Shape constructor
  // Args:
  //  worldline: worldline of a point in restframe of the shape
  //  o: a point in the shape, specified in rest frame
  //  a, b: a pair of vectors spanning the plane of the shape, specified in rest
  //   frame
  //  mat: material of the shape
  Shape2D(Line worldline, Vec3 o, Vec3 a, Vec3 b, Material mat = Material())
      : Object(worldline, mat), o_(o), a_(a), b_(b), normal_(NormalTo(a, b)) {}
  // Returns t such that x_0 + vel * t intersects the plane through the origin
  // with given normal if such t exists, and std::nullopt otherwise.
  static std::optional<float> IntersectionTime(Vec3 x_0, Vec3 vel, Vec3 plane_o,
                                               Vec3 plane_normal);

  // Given vectors x, o, a, b, returns u,v such that
  // x' = u * a + v * b is the closest approximation to x
  // in the plane going through o spanned by a and b
  // (mostly used when x is on the plane spanned by a and b)
  static UVCoordinates UVCoords(Vec3 x, Vec3 o, Vec3 a, Vec3 b);

  // Return ReferenceFrameHit for 2D shape, automating
  // inside/outside detection
  static ReferenceFrameHit MakeRfHit(Vec4 pos, Vec3 ray_vel,
                                     Vec3 outside_normal);

 private:
  // all in rest frame:
  const Vec3 o_;       // origin of shape in rest frame
  const Vec3 a_, b_;   // two vectors spanning the plane of the shape
  const Vec3 normal_;  // outside normal of the shape, (a ∧ b).Normalized()

  OptionalReferenceFrameHit intersect_in_rest_frame(const Line &) const;

  // Wrapper for is_inside_shape unpacking a UVCoordinates object into u and v.
  bool is_inside_shape(UVCoordinates uv) const {
    return is_inside_shape(uv.u, uv.v);
  }

  // Given UV coordinates u, v, specifying a point x = o + u*a + v*b in the
  // plane of the object, return true if the point x is inside the shape. To be
  // overridden by explicit shapes deriving from this class.
  virtual bool is_inside_shape(float u, float v) const = 0;
};

// Triangle specialisation of Shape2D class.
class Triangle : public Shape2D {
 public:
  // Triangle constructor
  // Args:
  //  worldline: worldline of coordinate origin
  //  A_rf, B_rf, C_rf: the three vertices of the triangle specified in
  //   restframe
  //  mat: triangle material
  Triangle(Line worldline, Vec3 A_rf, Vec3 B_rf, Vec3 C_rf,
           Material mat = Material())
      : Shape2D(worldline, A_rf, B_rf - A_rf, C_rf - A_rf, mat) {}
  // Triangle constructor
  // Args:
  //  O: 3-position of coordinate origin
  //  A_rf, B_rf, C_rf: triangle vertices in rest frame of O
  //  vel: velocity of O (in std frame)
  //  t: time at which O was given (in std frame)
  //  mat: triangle material
  Triangle(Vec3 O, Vec3 A_rf, Vec3 B_rf, Vec3 C_rf, Vec3 vel = kZero3,
           float t = 0, Material mat = Material())
      : Triangle(Line(Vec4(t, O), vel), A_rf, B_rf, C_rf, mat) {}

 private:
  // check if u > 0, v > 0, u+v <= 1
  bool is_inside_shape(float u, float v) const override;
};

// Parallelogram specialisation of Shape2D class.
class Parallelogram : public Shape2D {
 public:
  // Parallelogram constructor
  // Args:
  //  worldline: worldline of the coordinate origin
  //  o_rf: position of a vertex, given in rest frame
  //  a_rf, b_rf: vectors spanning the parallelogram, starting from o_rf. Given
  //   in rest frame.
  //  mat: parallelogram material
  Parallelogram(Line worldline, Vec3 o_rf, Vec3 a_rf, Vec3 b_rf,
                Material mat = Material())
      : Shape2D(worldline, o_rf, a_rf, b_rf, mat) {}
  // Parallelogram constructor
  // Args:
  //  worldline_c: worldline of center of parallelogram
  //  half_a, half_b: vectors from center to edge centers in rest frame
  //  mat: parallelogram material
  Parallelogram(Line worldline_c, Vec3 half_a, Vec3 half_b,
                Material mat = Material())
      : Parallelogram(worldline_c, -(half_a + half_b), 2 * half_a, 2 * half_b,
                      mat) {}
  // Parallelogram constructor
  // Args:
  //  C: position of center of parallelogram, given in standard frame
  //  half_a, half_b: vectors from center to edge centers in rest frame
  //  vel: velocity of C (in standard frame)
  //  t: time at which C is given (in standard frame)
  //  mat: parallelogram material
  Parallelogram(Vec3 C, Vec3 half_a, Vec3 half_b, Vec3 vel = kZero3,
                float t = 0.0f, Material mat = Material())
      : Parallelogram(Line(Vec4(t, C), vel), -(half_a + half_b), 2 * half_a,
                      2 * half_b, mat) {}
  // Parallelogram constructor (stationary)
  // Args:
  //  C: position of center of parallelogram, given in standard frame
  //  half_a, half_b: vectors from center to edge centers in rest frame
  //  mat: parallelogram material
  Parallelogram(Vec3 C, Vec3 half_a, Vec3 half_b, Material mat = Material())
      : Parallelogram(Line(Vec4(0, C), kZero3), -(half_a + half_b), 2 * half_a,
                      2 * half_b, mat) {}

 private:
  // check if u > 0, v > 0, u < 1, v < 1
  bool is_inside_shape(float u, float v) const override;
};

// Plane specialisation of Shape2D class. An infinite plane.
class Plane : public Shape2D {
 public:
  // Construct plane from worldline of a point
  // and two vectors `a_rf`, `b_rf` spanning the plane
  // in the restframe of that point
  // Plane constructor:
  // Args:
  //  worldline: worldline of a point on the plane
  //  a_rf, b_rf: vectors spanning the plane, in rest frame
  //  mat: plane material
  Plane(Line worldline, Vec3 a_rf, Vec3 b_rf, Material mat = Material())
      : Shape2D(worldline, kZero3, a_rf, b_rf, mat) {}
  // Plane constructor
  // Args:
  //  O: position of point on plane, given in standard frame
  //  a_rf, b_rf: vectors spanning plane, in rest frame
  //  vel: velocity of O (in standard frame)
  //  t: time at which O is given (in standard frame)
  //  mat: plane material
  Plane(Vec3 O, Vec3 a_rf, Vec3 b_rf, Vec3 vel = kZero3, float t = 0,
        Material mat = Material())
      : Plane(Line(Vec4(t, O), vel), a_rf, b_rf, mat) {}
  // Plane constructor (stationary)
  // Args:
  //  O: position of point on plane, given in standard frame
  //  a_rf, b_rf: vectors spanning the plane in rest frame
  //  mat: plane material
  Plane(Vec3 O, Vec3 a_rf, Vec3 b_rf, Material mat = Material())
      : Plane(O, a_rf, b_rf, kZero3, 0, mat) {}

 private:
  // Always returns true - the plane is infinite.
  bool is_inside_shape(float u, float v) const override;
};

// Box object, made up from 6 parallelogram.
// The six sides may be accessed individually, for example to
// change their materials.
class Box {
 public:
  // we say `up` & `down` to avoid the Bottom-Back intials clash
  Parallelogram right;
  Parallelogram back;
  Parallelogram up;
  Parallelogram left;
  Parallelogram front;
  Parallelogram down;

  // Box constructor
  // Args:
  //  worldline_c: worldline of point c in Box
  //  c_to_right, c_to_back, c_to_up: vector from c to center of right/back/top
  //   side, specified in rest frame
  //  mat: box material
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

  // Box constructor
  // Args:
  //  C: position of point c in box, in standard frame
  //  c_to_right, c_to_back, c_to_up: vector from c to center of right/back/top
  //   side, specified in rest frame
  //  vel: velocity of C
  //  t: time at which C was specified
  //  mat: box material
  Box(Vec3 C, Vec3 c_to_right, Vec3 c_to_back, Vec3 c_to_up, Vec3 vel = kZero3,
      float t = 0.0f, Material mat = Material())
      : Box(Line(Vec4(t, C), vel), c_to_right, c_to_back, c_to_up, mat) {}

  // Box constructor (stationary)
  // Args:
  //  C: position of point c in box, in standard frame
  //  c_to_right, c_to_back, c_to_up: vector from c to center of right/back/top
  //   side, specified in rest frame
  //  mat: box material
  Box(Vec3 C, Vec3 c_to_right, Vec3 c_to_back, Vec3 c_to_up, Material mat)
      : Box(C, c_to_right, c_to_back, c_to_up, kZero3, 0, mat) {}
};

// Camera class - at least one must be added to each Scene.
class Camera {
 public:
  // Constructs camera from `worldline` of camera sensor and screen data (given
  // in the restframe of the camera with the camera at the origin).
  // Pixels are assumed to be square in the camera frame.
  // Args:
  //  worldline: worldline of (point-like) camera sensor
  //  cam_to_screen_center: vector from sensor to screen center in camera rest
  //   frame
  //  screen_center_to_right: vector from screen center to right edge of screen
  //   in camera rest frame
  //  screen_width: screen width in pixels
  //  screen_height: screen height in pixels
  Camera(Line worldline, Vec3 cam_to_screen_center, Vec3 screen_center_to_right,
         int width, int height)
      : worldline_(worldline),
        screen_center_(cam_to_screen_center),
        screen_center_to_right_(screen_center_to_right),
        screen_width_(width),
        screen_height_(height) {}

  // Constructs camera from `worldline` of camera sensor and screen data (given
  // in the restframe of the camera with the camera at the origin).
  // Pixels are assumed to be square in the camera frame.
  // Args:
  //  worldline: worldline of (point-like) camera sensor
  //  cam_to_screen_center: vector from sensor to screen center in camera rest
  //   frame
  //  screen_center_to_right: vector from screen center to right edge of screen
  //   in camera rest frame
  //  screen_width: screen width in pixels
  //  screen_height: screen height in pixels
  Camera(Vec3 origin, Vec3 cam_to_screen_center, Vec3 screen_center_to_right,
         int width, int height, Vec3 vel = kZero3, float t = 0)
      : Camera(Line(Vec4(t, origin), vel), cam_to_screen_center,
               screen_center_to_right, width, height) {}

  int GetScreenWidth() const { return screen_width_; }
  int GetScreenHeight() const { return screen_height_; }
  Vec3 GetVelocity() const { return worldline_.vel; }
  // Get light rays that make up the image at specified time `t_cam_frame`,
  // given in the restframe of the camera.
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

// Scene class. Contains the data for a scene to be rendered: a number of
// uniformly moving objects and a list of cameras.
// Note: cameras are not rendered.
class Scene {
 public:
  // Polymorphic method to add an object or a camera to the Scene.
  template <class Derived_Object>
  void Add(const Derived_Object &_obj) {
    objects_.push_back(std::make_shared<Derived_Object>(_obj));
  }

  // Specialisation to Box: this is implemented as adding the six sides.
  void Add(const Box &box) {
    Add(box.front);
    Add(box.down);
    Add(box.left);
    Add(box.back);
    Add(box.up);
    Add(box.right);
  }

  // Specialisation to Camera: adds `cam` to camera list.
  void Add(const Camera &cam) { cameras_.push_back(cam); }

  // Get camera with given index.
  const Camera &GetCamera(int index = 0) const;

  // Finds the nearest (most recent) intersection of ray with any object
  // in the scene.
  // Returns HitRecord if an intersection exists and std::nullopt otherwise.
  OptionalHitRecord MostRecentHit(const Line &ray) const;

  // Color spectrum of a ray which came from the background (without hitting an
  // object in the scene).
  Spectrum EscapedRayColor(const Line &ray) const;

  // Set color spectrum of diffuse ambient background light.
  Scene &SetAmbientBackground(const Spectrum &bg_spectrum);

  // Set a directional ambient light source, at rest in standard frame.
  // Args:
  //  skylight: color spectrum of light
  //  skylight_direction: direction of light coming from skylight
  //  skylight_exponent: exponent of falloff. When this is equal to 1,
  //  brightness falls off with the cosine law.
  Scene &SetSkylight(const Spectrum skylight, Vec3 skylight_direction,
                     float skylight_exponent = 1.0f);

 private:
  CameraList cameras_;
  ObjectList objects_;
  Spectrum ambient_background_ = Spectrum::Black;
  Spectrum skylight_ = Spectrum::Black;
  Vec3 skylight_direction_ = kDefaultNormal;
  float skylight_exponent_ = 1.0f;
};

#endif
