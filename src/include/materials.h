// Copyright 2023 Christoph Weis
#ifndef MATERIALS_H_INCLUDED
#define MATERIALS_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "colors.h"
#include "math.h"

const unsigned int kStandardSeed = 42;

// Material class, implemented non-polymorphically
class Material {
 public:
  // Material
  // Args:
  //  albedo: overall brightness multiplier applied to all scattered light
  //  diffuse: fraction of light to be scattered diffusely (Lambertian model)
  //  specular: fraction of light to be reflected (for Mirrors & highlights)
  //  fuzz: disturbance applied to reflected light (length of random vector
  //   added to reflected velocity)
  //  absorption: color Spectrum absorbed, treated as a fraction (!)
  //  emission: color Spectrum emitted, treated as specific intensity
  //  dielectric: bool controlling if medium is dielectric
  //   if true: all parameters but dielectric_n and emission are ignored
  //  dielectric_n: dielectric index of medium, ratio of vacuum light speed
  //   to light speed in dielectric medium (huge simplification: the index is
  //   assumed constant across wavelengths)
  //  subsurface_scatters: fraction of light to be scattered by the
  //   subsurface (implemented as simply a second Lambertian scatter)
  //  subsurface_absorption: color absorption spectrum of subsurface layer
  //   (treated as fraction)
  Material(float albedo = 1.0f, float diffuse = 1.0f, float specular = 0.0f,
           float fuzz = 0.0f, Spectrum absorption = Spectrum::Black,
           Spectrum emission = Spectrum::Black, bool dielectric = false,
           float dielectric_n = 1.0f, float subsurface_scatters = 0.0f,
           Spectrum subsurface_absorption = Spectrum::Black)
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

  // Useful Material Types
  static Material Dielectric(float dielectric_n);
  static Material Metal(float albedo = 1.0f,
                        Spectrum absorption = Spectrum::Black,
                        float fuzz = 0.0f);
  static Material Mirror(float fuzz = 0.0f);
  static Material Light(Spectrum emission = Spectrum::White,
                        Spectrum absorption = Spectrum::White);

  static const Material Water;
  static const Material Glass;
  static const Material Diamond;

  // Returns emission color spectrum (treated as specific intensity)
  Spectrum EmissionSpectrum(const ReferenceFrameHit &rf_hit) const;

  // Returns color absorption spectrum (treated as fraction)
  Spectrum AbsorptionCurve(const ReferenceFrameHit &rf_hit) const;

  // Performs inverse scatter against material given reference frame hit data.
  // Returns velocity of scattered ray, changes cumulative transform to reflect
  // transformation of light in brightness and frequency. This incorporates
  // both absorption effects and special relativistic transformation.
  Vec3 InverseScatter(const ReferenceFrameHit &rf_hit,
                      SpectrumTransform *cumulative_transform) const;

  // Property setter for float properties
  // Allowed names:
  // "albedo", "diffuse", "specular", "fuzz", "dielectric_n",
  // "subsurface_scatters"
  Material &SetProperty(std::string scalar_property_name, float new_value);

  // Property setter for Spectrum properties
  // Allowed names:
  // "emission", "absorption", "subsurface_absorption"
  Material &SetProperty(std::string spectral_property_name, Spectrum new_value);

  // Property setter for bool properties
  // Allowed names:
  // "dielectric"
  Material &SetProperty(std::string bool_property_name, bool new_value);

 private:
  float albedo_ = 1.0f;        // brightness multiplier
  float diffuse_ = 1.0f;       // lambertian scatter multiplier
  float specular_ = 0.0f;      // specular reflection multiplier
  float fuzz_ = 0.0f;          // fuzz modifier (jitters specular reflection)
  bool dielectric_ = false;    // if medium is dielectric, only n is used
  float dielectric_n_ = 1.0f;  // refractive index
  float subsurface_scatters_ = 0.0f;  // subsurface multiplier - assumed diffuse
  Spectrum emission_ = Spectrum::Black;  // spectrum emitted - always applied
  Spectrum absorption_ = Spectrum::Black;  // spectrum absorbed - diffuse & specular
  Spectrum subsurface_absorption_ =
      Spectrum::Black;  // spectrum absorbed - subsurface

  float non_dielectric_multipliers_() const {
    return diffuse_ + specular_ + subsurface_scatters_;
  }

  Vec3 DielectricInverseScatter(const ReferenceFrameHit &rf_hit,
                                SpectrumTransform *cumulative_transform) const;
  Vec3 LambertianInverseScatter(const ReferenceFrameHit &rf_hit,
                                SpectrumTransform *cumulative_transform) const;
  Vec3 SpecularInverseScatter(const ReferenceFrameHit &rf_hit,
                              SpectrumTransform *cumulative_transform) const;
  Vec3 SubsurfaceInverseScatter(const ReferenceFrameHit &rf_hit,
                                SpectrumTransform *cumulative_transform) const;
};

#endif
