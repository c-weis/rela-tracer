// Copyright 2023 Christoph Weis
#include "include/scene.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

#include "include/colors.h"
#include "include/materials.h"
#include "scene.h"

// Constant allowing for rounding errors in
// computing intersections of scattered rays
const float kTimeEpsilon = 5e-6f;

/*
  CAMERA DATA
*/
// Computes the list of image rays for the image received at
// the camera after delta_t has elapsed
// (counted in camera ref frame, from camera origin time)
LineList Camera::ImageRays(float delta_t) const {
  float aspect = static_cast<float>(screen_width_) / screen_height_;
  Vec3 central_line = screen_centre_;
  Vec3 right = screen_right_;
  Vec3 down = Cross(central_line, right) / aspect;

  LineList rays;
  rays.reserve(screen_width_ * screen_height_);
  const float dx = 2.0f / screen_width_;
  const float dy = 2.0f / screen_height_;
  for (int j = 0; j < screen_height_; j++) {
    for (int i = 0; i < screen_width_; i++) {
      Vec2 random_vec = RandomVectorInUnitDisk();
      float x = -1.0f + (i + 0.5f + random_vec.x * 0.7071) * dx;
      float y = -1.0f + (j + 0.5f + random_vec.y * 0.7071) * dy;
      Line ray_camera_frame =
          Line(worldline_.PosAfter(delta_t),
               -(central_line + x * right + y * down).NormalizedNonzero());
      // transform to standard frame
      rays.push_back(ray_camera_frame.TransformedFromFrame(worldline_));
    }
  }

  return rays;
}

Scene &Scene::AddCamera(const Camera &camera) {
  cameras_.push_back(camera);
  return *this;
}

const Camera &Scene::GetCamera(int index) const {
  if (index > cameras_.size() - 1) {
    throw std::invalid_argument("Camera index " + std::to_string(index) +
                                " exceeds number of cameras in Scene.");
  }
  return cameras_.at(index);
}

Scene &Scene::SetSkylight(Spectrum skylight_spectrum, Vec3 skylight_direction) {
  skylight_ = skylight_spectrum;
  skylight_direction_ = skylight_direction;
  return *this;
}

/*
  HITRECORD THINGS
*/

inline float HitRecord::hit_time() const {
  return (rf.pos + obj->GetOrigin()).TransformedFromFrame(obj->GetVelocity()).t;
}

/*
  ABSTRACT OBJECT MOVEMENT AND INTERSECTION
*/
Vec4 Object::PosAt(float time) const { return worldline_.PosAt(time); }

Vec4 Object::PosAfter(float d_time) const {
  return worldline_.PosAfter(d_time);
}

OptionalHitRecord Object::intersect(const Line &ray) const {
  OptionalReferenceFrameHit rf_hit =
      intersect_in_rest_frame(ray.TransformedToFrame(worldline_));
  if (!rf_hit.has_value()) {
    return std::nullopt;
  }
  return HitRecord(rf_hit.value(), this);
}

std::optional<float> plane_intersection_time(Vec3 x_0, Vec3 vel,
                                             Vec3 plane_normal) {
  // assume plane goes through the origin
  // TODO(c): build in division safety
  // alternatively return no hit here
  float v_n = Dot3(vel, plane_normal);
  if (v_n >= 0 && v_n < kDivisionEpsilon) {
    // v_n = kDivisionEpsilon;
    return std::nullopt;
  } else if (v_n < 0 && v_n > -kDivisionEpsilon) {
    // v_n = -kDivisionEpsilon;
    return std::nullopt;
  }
  float delta = -Dot3(x_0, plane_normal) / v_n;
  // account for fact that rays travel BACKWARD in time
  if (delta >= -kTimeEpsilon) return std::nullopt;
  return delta;
}

UVCoordinates uv_coordinates(Vec3 x, Vec3 a, Vec3 b) {
  float a_b = Dot3(a, b);
  float a_a = a.NormSq();
  float b_b = b.NormSq();
  float x_a = Dot3(x, a);
  float x_b = Dot3(x, b);
  float det = a_a * b_b - a_b * a_b;
  if (det < kDivisionEpsilon) {
    // a and b are collinear - throw exception
    throw std::invalid_argument(
        "The given vectors a and b are collinear, "
        "the plane spanned by them is degenerate.");
  }

  return UVCoordinates{.u = (b_b * x_a - a_b * x_b) / det,
                       .v = (-a_b * x_a + a_a * x_b) / det};
}

ReferenceFrameHit rf_hit_from(Vec4 pos, Vec3 ray_vel, Vec3 outside_normal) {
  // Checks if ray_vel is aligned or anti-aligned with outside_normal,
  // deduces whether we hit in- or outside
  if (Dot3(ray_vel, outside_normal) < 0) {
    // we hit the inside: invert normal, specify outside=false
    return ReferenceFrameHit(pos, -outside_normal, ray_vel, false);
  }
  return ReferenceFrameHit(pos, outside_normal, ray_vel);
}

/*
  SPECIFIC OBJECTS
*/

float Sphere::rad() { return rad_; }

OptionalReferenceFrameHit Sphere::intersect_in_rest_frame(
    const Line &ray) const {
  Vec3 x_0 = ray.origin.r;
  Vec3 v = ray.vel;

  float proj = Dot3(v, x_0);
  float v_sq = v.NormSq();

  // check if there are any intersections, return empty if not
  float discr = proj * proj - v_sq * (x_0.NormSq() - rad_ * rad_);
  if (discr < 0) return std::nullopt;

  // calculate the two intersection times
  float delta_minus = (-proj - sqrtf(discr)) / v_sq;
  float delta_plus = (-proj + sqrtf(discr)) / v_sq;

  // account for fact that rays travel BACKWARD in time
  // no intersection if 0 <= delta_minus <= delta_plus
  if (delta_minus >= -kTimeEpsilon) return std::nullopt;

  if (delta_plus >= -kTimeEpsilon) {
    // delta_minus <= 0 <= delta_plus
    // we hit the INSIDE of the sphere at delta_minus
    Vec4 pos = ray.PosAfter(delta_minus);

    return ReferenceFrameHit(pos, -pos.r / rad_, v, false);
  }
  // delta_minus <= delta_plus <= 0
  // we hit the OUTSIDE of the sphere at delta_plus
  Vec4 pos = ray.PosAfter(delta_plus);

  return ReferenceFrameHit(pos, pos.r / rad_, v);
}

// Finds nearest hit of the ray with an object in the scene.
// We compute the `hit_time` as experienced in the standard
// inertial frame and return the hit corresponding to the
// earliest time.
OptionalHitRecord Scene::MostRecentHit(const Line &ray) const {
  OptionalHitRecord most_recent_hit;
  // positive `best time` is code for no hit
  // rays travel BACKWARD in time
  float best_time = 1;

  // Loop over objects, find nearest (=most recent) intersection
  for (const auto &obj : objects_) {
    OptionalHitRecord obj_hit = obj->intersect(ray);
    if (!obj_hit) continue;
    float obj_time = obj_hit.value().hit_time();
    if (obj_time > -kTimeEpsilon) continue;

    if (best_time > 0 || (obj_time > best_time)) {
      most_recent_hit = obj_hit;
      best_time = obj_time;
    }
  }

  return most_recent_hit;
}

Spectrum Scene::EscapedRayColor(const Line &ray) const {
  return ambient_background_ + skylight_ * Dot3(skylight_direction_, ray.vel);
}

Scene &Scene::SetAmbientBackground(const Spectrum &bg_spectrum) {
  ambient_background_ = bg_spectrum;
  return *this;
}

OptionalReferenceFrameHit Shape2D::intersect_in_rest_frame(
    const Line &ray) const {
  Vec3 x_0 = ray.origin.r;
  Vec3 vel = ray.vel;

  std::optional<float> delta_ = plane_intersection_time(x_0, vel, normal_);
  if (!delta_.has_value()) return std::nullopt;
  float delta = delta_.value();

  Vec4 pos = ray.PosAfter(delta);
  UVCoordinates uv = uv_coordinates(pos.r, a_, b_);

  if (!is_inside_shape(uv)) return std::nullopt;

  return rf_hit_from(pos, vel, normal_);
}

bool Plane::is_inside_shape(float u, float v) const { return true; }

bool Triangle::is_inside_shape(float u, float v) const {
  if (u < 0 || v < 0 || u + v > 1) {
    return false;
  }
  return true;
}

bool Parallelogram::is_inside_shape(float u, float v) const {
  if (u < 0 || u > 1 || v < 0 || v > 1) {
    return false;
  }
  return true;
}
