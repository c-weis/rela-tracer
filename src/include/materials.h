// Copyright 2023 Christoph Weis
#ifndef MATERIALS_H_INCLUDED
#define MATERIALS_H_INCLUDED

#include <algorithm>
#include <memory>
#include <random>
#include <vector>

#include "colors.h"
#include "math.h"

const uint kStandardSeed = 42;

struct ScatterData {
  // struct for a weighted light ray,
  // plus an auxiliary material flag
  // to be used for computation of
  // scattered color
  ScatterData(Vec3 vel_, float weight_ = 1.0f, int material_data_ = 0)
      : vel(vel_), weight(weight_), material_data(material_data_) {}
  Vec3 vel;
  float weight;
  int material_data;
};

inline ScatterData operator*(const ScatterData &s_data, float factor) {
  return ScatterData(s_data.vel, s_data.weight * factor, s_data.material_data);
}

class Component {
 public:
  virtual std::vector<ScatterData> InverseScatter(
      const ReferenceFrameHit &rf_hit, int scatter_rays,
      std::mt19937 *random_number_generator) const {
    return {};
  }

  virtual ColorData ScatteredColor(const ColorData &pre_scatter_color,
                                   const ReferenceFrameHit &rf_hit,
                                   const Vec3 &pre_scattered,
                                   int material_data = 0) const {
    return kBlack;
  }

  virtual ColorData EmittedColor(const ReferenceFrameHit &rf_hit) const {
    return kBlack;
  }
};

struct WeightedComponent {
  template <class ComponentType>
  WeightedComponent(const std::shared_ptr<ComponentType> &comp_,
                    float weight_ = 1.0f)
      : comp(comp_), weight(weight_) {}

  template <class ComponentType>
  WeightedComponent(const ComponentType &comp_, float weight_ = 1.0f)
      : WeightedComponent(std::make_shared<ComponentType>(comp_), weight_) {}

  std::shared_ptr<Component> comp;
  float weight;
};

/*
inline WeightedComponent operator*(const WeightedComponent &w_comp,
                                   float factor) {
  return WeightedComponent(w_comp.comp, w_comp.weight * factor);
}

inline WeightedComponent operator*(float factor,
                                   const WeightedComponent &w_comp) {
  return factor * w_comp;
}
*/

typedef std::vector<WeightedComponent> ComponentList;
class Material : public Component {
 public:
  Material(WeightedComponent w_comp) : components_({w_comp}) {}

  template <class ComponentType>
  Material(const ComponentType &comp) : Material(WeightedComponent(comp)) {}

  template <class... ComponentTypes>
  Material(WeightedComponent w_comp, ComponentTypes... more_comps)
      : Material(w_comp) {
    Add(more_comps...);
  }

  template <class ComponentType, class... ComponentTypes>
  Material(ComponentType comp1, ComponentTypes... more_comps)
      : Material(WeightedComponent(comp1)) {
    Add(more_comps...);
  }

  void Add(WeightedComponent w_comp) { components_.push_back(w_comp); }

  template <class ComponentType>
  void Add(ComponentType comp) {
    components_.push_back(WeightedComponent(comp));
  }

  template <class ComponentType, class... ComponentTypes>
  void Add(ComponentType comp1, ComponentTypes... more_comps) {
    Add(comp1);
    Add(more_comps...);
  }

  friend Material operator+(const Material &mat1, const Material &mat2) {
    Material new_mat(mat1);
    new_mat.components_.insert(new_mat.components_.end(),
                               mat2.components_.begin(),
                               mat2.components_.end());
    return new_mat;
  }

  friend Material operator*(const Material &mat, float factor) {
    Material new_mat(mat);
    std::for_each(
        new_mat.components_.begin(), new_mat.components_.end(),
        [factor](WeightedComponent &w_comp) { w_comp.weight *= factor; });
    return new_mat;
  }

  std::vector<ScatterData> InverseScatter(
      const ReferenceFrameHit &rf_hit, int scatter_rays,
      std::mt19937 *random_number_generator) const override;

  ColorData ScatteredColor(const ColorData &pre_scatter_color,
                           const ReferenceFrameHit &rf_hit,
                           const Vec3 &pre_scattered,
                           int material_data) const override;

  ColorData EmittedColor(const ReferenceFrameHit &rf_hit) const override;

 private:
  ComponentList components_;
};

inline Material operator+(const WeightedComponent &comp1,
                          const WeightedComponent &comp2) {
  return Material(comp1, comp2);
}

inline Material operator*(float factor, const Material &mat) { return mat * factor; }

inline Material operator+(const Material &mat, const WeightedComponent &comp) {
  Material new_comp_mat(mat);
  new_comp_mat.Add(comp);
  return new_comp_mat;
}

class Lambertian : public Component {
 public:
  Lambertian(float albedo = 1.0f) : albedo_(albedo) {}

  std::vector<ScatterData> InverseScatter(
      const ReferenceFrameHit &rf_hit, int scatter_rays,
      std::mt19937 *random_number_generator) const override;

  ColorData ScatteredColor(const ColorData &pre_scatter_color,
                           const ReferenceFrameHit &rf_hit,
                           const Vec3 &pre_scattered,
                           int material_data) const override;

 private:
  float albedo_;
};

class MonochromaticLight : public Component {
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

  ColorData EmittedColor(const ReferenceFrameHit &rf_hit) const override;

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

class Dielectric : public Component {
 public:
  Dielectric(float n) : n_(n) {}

  std::vector<ScatterData> InverseScatter(
      const ReferenceFrameHit &rf_hit, int scatter_rays,
      std::mt19937 *random_number_generator) const override;

  ColorData ScatteredColor(const ColorData &pre_scatter_color,
                           const ReferenceFrameHit &rf_hit,
                           const Vec3 &pre_scattered,
                           int material_data) const override;

 private:
  float n_;  // refractive index
};

const Dielectric kWater(1.3f);
const Dielectric kGlass(1.5f);
const Dielectric kDiamond(2.4f);

class Metal : public Component {
 public:
  // TODO(c): build in clamps for the values
  Metal(float albedo = 1.0f, float fuzz = 0.0f)
      : albedo_(albedo), fuzz_(fuzz) {}

  std::vector<ScatterData> InverseScatter(
      const ReferenceFrameHit &rf_hit, int scatter_rays,
      std::mt19937 *random_number_generator) const override;

  ColorData ScatteredColor(const ColorData &pre_scatter_color,
                           const ReferenceFrameHit &rf_hit,
                           const Vec3 &pre_scattered,
                           int material_data) const override;

 private:
  float albedo_;
  float fuzz_;
};

const Metal kMirror(1.0f, 0.0f);

#endif