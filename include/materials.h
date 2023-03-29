// Copyright 2023 Christoph Weis
#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include <vector>

#include "include/math.h"

const uint kStandardSeed = 42;

/*
  COLOR STRUCTURES
*/
struct RGBData {
  int R, G, B;
  static RGBData gray(int shade);
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
  LightRay TransformedToFrame(const Vec3 &) const;
  LightRay TransformedFromFrame(const Vec3 &) const;
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

  rgbData To_rgb() const;
  RGBData ToRGB() const;
};

// Method culled because it is inefficient
// ColorData operator+(const ColorData &, const ColorData &);
ColorData operator*(const ColorData &, float);
ColorData operator*(float, const ColorData &);
ColorData operator/(const ColorData &, float);

const ColorData kBlack;
// const ColorData kWhite(1.0f);

/*
  MATERIALS SECTION
*/
class Scatterer {
 public:
  virtual std::vector<Vec3> InverseScatter(
      Vec3 scattered_ray, Vec3 object_normal, int max_number_of_rays,
      std::mt19937 *random_number_generator) const = 0;

  virtual const ColorData ScatteredColor(
      const ReferenceFrameHit &rf_hit, const Vec3 &pre_scattered,
      const ColorData &pre_scatter_color) const = 0;

 private:
  ColorData color_;
};

class Emitter {
 public:
  virtual const ColorData EmittedColor(
      const ReferenceFrameHit &rf_hit) const = 0;
};

class Material {
 public:
  template <class ScatterType, class EmitterType>
  Material(const ScatterType *scatterer, const EmitterType *emitter)
      : scatterer_(scatterer), emitter_(emitter) {}
  template <class EmitterType>
  Material(std::nullptr_t, const EmitterType *emitter)
      : scatterer_(nullptr), emitter_(emitter) {}
  template <class ScatterType>
  Material(const ScatterType *scatterer, std::nullptr_t)
      : scatterer_(scatterer), emitter_(nullptr) {}

  Material &SetScatterer(Scatterer *scatterer);
  Material &SetEmitter(Emitter *emitter);

  // call scatterer
  std::vector<Vec3> InverseScatter(Vec3 scattered_ray, Vec3 object_normal,
                                   int max_number_of_rays,
                                   std::mt19937 *random_number_generator) const;

  // call scatterer
  const ColorData ScatteredColor(const ReferenceFrameHit &rf_hit,
                                 const Vec3 &pre_scattered,
                                 const ColorData &pre_scatter_color) const;

  // call emitter
  const ColorData EmittedColor(const ReferenceFrameHit &rf_hit) const;

 private:
  const Scatterer *scatterer_;
  const Emitter *emitter_;
};

class LambertianScatter : public Scatterer {
 public:
  std::vector<Vec3> InverseScatter(
      Vec3 scattered_ray, Vec3 object_normal, int max_number_of_rays,
      std::mt19937 *random_number_generator) const override;

  const ColorData ScatteredColor(
      const ReferenceFrameHit &rf_hit, const Vec3 &pre_scattered,
      const ColorData &pre_scatter_color) const override;

  static Vec3 RandomUnitVector(std::mt19937 *r_gen);
  static Vec3 LambertianInverseRay(Vec3 object_normal, std::mt19937 *r_gen);

 private:
  /*
    TODO(c): define absorption data, use in Scattered Color
  */
};

class MonochromaticLight : public Emitter {
 public:
  /*
    Specify monochromatic light source by giving brightness
    (`brightness`) and angular frequency `omega`. The corresponding
    wavenumber in a medium M is k_M = omega / c_M, where
    c_M is the speed of light in M.
    Unlike the wavenumber, `omega` is constant across media interfaces.
    This way of specifying the light source allows us to place light sources
    inside non-vacuum media.

    We use units of nanometers for wavelengths and set c_vac=1, which forces us
    to use units of `radians/nanosecond` for `omega`.
  */
  MonochromaticLight(float omega, float brightness = 1.0f)
      : brightness_(brightness), omega_(omega) {}
  /*
    Calculates angular frequency from (vacuum) wavelength of given light.
  */
  static MonochromaticLight FromWavelength(float lambda_vacuum,
                                           float brightness = 1.0f) {
    return MonochromaticLight(2 * M_PI / lambda_vacuum, brightness);
  };

  const ColorData EmittedColor(const ReferenceFrameHit &rf_hit) const override;

 private:
  const float brightness_;  // think: photons per second
  const float omega_;       // angular frequency
};
const MonochromaticLight kBrightRed =
    MonochromaticLight::FromWavelength(680, 1.0f);
const MonochromaticLight kBrightOrange =
    MonochromaticLight::FromWavelength(605, 1.0f);
const MonochromaticLight kBrightYellow =
    MonochromaticLight::FromWavelength(580, 1.0f);
const MonochromaticLight kBrightGreen =
    MonochromaticLight::FromWavelength(530, 1.0f);
const MonochromaticLight kBrightCyan =
    MonochromaticLight::FromWavelength(490, 1.0f);
const MonochromaticLight kBrightBlue =
    MonochromaticLight::FromWavelength(460, 1.0f);
const MonochromaticLight kBrightPurple =
    MonochromaticLight::FromWavelength(410, 1.0f);
const MonochromaticLight kLightRed =
    MonochromaticLight::FromWavelength(680, 0.5f);
const MonochromaticLight kLightOrange =
    MonochromaticLight::FromWavelength(605, 0.5f);
const MonochromaticLight kLightYellow =
    MonochromaticLight::FromWavelength(580, 0.5f);
const MonochromaticLight kLightGreen =
    MonochromaticLight::FromWavelength(530, 0.5f);
const MonochromaticLight kLightCyan =
    MonochromaticLight::FromWavelength(490, 0.5f);
const MonochromaticLight kLightBlue =
    MonochromaticLight::FromWavelength(460, 0.5f);
const MonochromaticLight kLightPurple =
    MonochromaticLight::FromWavelength(410, 0.5f);
const MonochromaticLight kDimRed =
    MonochromaticLight::FromWavelength(680, 0.1f);
const MonochromaticLight kDimOrange =
    MonochromaticLight::FromWavelength(605, 0.1f);
const MonochromaticLight kDimYellow =
    MonochromaticLight::FromWavelength(580, 0.1f);
const MonochromaticLight kDimGreen =
    MonochromaticLight::FromWavelength(530, 0.1f);
const MonochromaticLight kDimCyan =
    MonochromaticLight::FromWavelength(490, 0.1f);
const MonochromaticLight kDimBlue =
    MonochromaticLight::FromWavelength(460, 0.1f);
const MonochromaticLight kDimPurple =
    MonochromaticLight::FromWavelength(410, 0.1f);

class Mirror : public Scatterer {
 public:
  std::vector<Vec3> InverseScatter(
      Vec3 scattered_ray, Vec3 object_normal, int max_number_of_rays,
      std::mt19937 *random_number_generator) const override;

  const ColorData ScatteredColor(
      const ReferenceFrameHit &rf_hit, const Vec3 &pre_scattered,
      const ColorData &pre_scatter_color) const override;

 private:
};

// TODO(c): add Material based on mirror
// class Metal : public Material;

// TODO(c): add refractive materials
// class Dielectric : public Material{
//};
