// Copyright 2023 Christoph Weis
#ifndef MATERIALS_H_INCLUDED
#define MATERIALS_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "colors.h"
#include "math.h"

const unsigned int kStandardSeed = 42;

struct ScatterData {
  // class for a light ray
  ScatterData(Vec3 vel) : vel_(vel) {}

  // default constructor
  ScatterData() : ScatterData(kZero3) {}

  Vec3 vel_;  // ray velocity
};

class Material {
 public:
  Material(float albedo = 1.0f, float diffuse = 1.0f, float specular = 0.0f,
           float fuzz = 0.0f, Spectrum absorption = kBlack,
           Spectrum emission = kBlack, bool dielectric = false,
           float dielectric_n = 1.0f, float subsurface_scatters = 0.0f,
           Spectrum subsurface_absorption = kBlack)
      : dielectric_(dielectric),
        dielectric_n_(dielectric_n),
        albedo_(albedo),
        diffuse_(diffuse),
        specular_(specular),
        fuzz_(fuzz),
        subsurface_scatters_(subsurface_scatters),
        emission_(emission),
        absorption_(absorption),
        subsurface_absorption_(subsurface_absorption) {}

  Spectrum EmissionSpectrum(const ReferenceFrameHit &rf_hit) const;

  Spectrum AbsorptionCurve(const ReferenceFrameHit &rf_hit) const;

  ScatterData InverseScatter(const ReferenceFrameHit &rf_hit,
                             SpectrumTransform &cumulative_transform) const;

  Material &SetProperty(std::string scalar_property_name, float new_value);
  Material &SetProperty(std::string spectral_property_name, Spectrum new_value);
  Material &SetProperty(std::string bool_property_name, bool new_value);

 private:
  float albedo_ = 1.0f;        // brightness multiplier
  float diffuse_ = 1.0f;       // lambertian scatter multiplier
  float specular_ = 0.0f;      // specular reflection multiplier
  float fuzz_ = 0.0f;          // fuzz modifier (jitters specular reflection)
  bool dielectric_ = false;    // if medium is dielectric, only n is used
  float dielectric_n_ = 1.0f;  // refractive index
  float subsurface_scatters_ = 0.0f;  // subsurface multiplier - assumed diffuse
  Spectrum emission_ = kBlack;        // spectrum emitted - always applied
  Spectrum absorption_ = kBlack;  // spectrum absorbed - diffuse and specular
  Spectrum subsurface_absorption_ = kBlack;  // spectrum absorbed - subsurface

  float non_dielectric_multipliers_() const {
    return diffuse_ + specular_ + subsurface_scatters_;
  }

  ScatterData DielectricInverseScatter(
      const ReferenceFrameHit &rf_hit,
      SpectrumTransform &cumulative_transform) const;
  ScatterData LambertianInverseScatter(
      const ReferenceFrameHit &rf_hit,
      SpectrumTransform &cumulative_transform) const;
  ScatterData SpecularInverseScatter(
      const ReferenceFrameHit &rf_hit,
      SpectrumTransform &cumulative_transform) const;
  ScatterData SubsurfaceInverseScatter(
      const ReferenceFrameHit &rf_hit,
      SpectrumTransform &cumulative_transform) const;
};

Material DielectricMaterial(float dielectric_n);
Material Metal(float albedo = 1.0f, Spectrum absorption = kBlack,
               float fuzz = 0.0f);
Material Mirror(float fuzz = 0.0f);

/*
  Specify monochromatic light source by giving brightness
  (`brightness`) and vacuum wavelength `lambda`. The corresponding
  angular frequency is 2 pi/lambda.
  Unlike the wavelength, `omega` is constant across media interfaces.

  We use units of nanometers for wavelengths and set c_vac=1, which forces us
  to use units of `radians/nanosecond` for `omega`.
*/

const Material kWater = DielectricMaterial(1.0f);
const Material kGlass = DielectricMaterial(1.5f);
const Material kDiamond = DielectricMaterial(2.4f);

#endif