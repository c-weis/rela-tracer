// Copyright 2023 Christoph Weis
#ifndef COLORS_H_INCLUDED
#define COLORS_H_INCLUDED

#include <vector>

#include "math.h"

/*
  COLOR STRUCTURES
*/
struct RGBData {
  int R, G, B;
};

struct rgbData {
  float r, g, b;
  RGBData ToRGB(float rescale_factor = 255.0f, int max = 255) const;
};
rgbData operator*(const rgbData &, float);
rgbData operator*(float, const rgbData &);
rgbData operator/(const rgbData &, float);

struct XYZData {
  float X, Y, Z;
  // compute XYZ from wavelength
  XYZData(float _X, float _Y, float _Z) : X(_X), Y(_Y), Z(_Z) {}
  XYZData(float lambda);

  rgbData To_rgb() const;
};

XYZData operator+(const XYZData &, const XYZData &);
XYZData operator*(const XYZData &, float);
XYZData operator*(float, const XYZData &);
XYZData operator/(const XYZData &, float);

struct LightRay {
  float brightness;
  Vec4 k;  // wave 4-vector (\omega, kx, ky, kz)

  LightRay(float photons_, Vec4 k_) : brightness(photons_), k(k_) {}
  LightRay &operator*=(float);
  LightRay &operator/=(float);
  // The wave vector transformed like any 4-vector,
  // while `brightness` should transform like 1/time.
  // (Geometric effects are already accounted for by
  // the statistical model we use to detect incoming light,
  // but the number of photons / second emitted by a light
  // source is dependent on the observer frame.)
  LightRay &TransformToFrame(const Vec3 &);
  LightRay &TransformFromFrame(const Vec3 &);
};
LightRay operator*(const LightRay &, float);
LightRay operator*(float, const LightRay &);

struct ColorData {
  std::vector<LightRay> light_rays;

  explicit ColorData(std::vector<LightRay> _light_rays = {})
      : light_rays(_light_rays) {}

  ColorData &operator+=(const ColorData &);
  ColorData &operator*=(float);
  ColorData &operator/=(float);

  ColorData &TransformToFrame(const Vec3 &);
  ColorData &TransformFromFrame(const Vec3 &);
  ColorData TransformedToFrame(const Vec3 &);
  ColorData TransformedFromFrame(const Vec3 &);

  rgbData To_rgb() const;
  RGBData ToRGB() const;
};

// Method culled because it is inefficient
// ColorData operator+(const ColorData &, const ColorData &);
ColorData operator*(const ColorData &, float);
ColorData operator*(float, const ColorData &);
ColorData operator/(const ColorData &, float);

static const ColorData kBlack;
// static const ColorData kWhite(1.0f);

#endif