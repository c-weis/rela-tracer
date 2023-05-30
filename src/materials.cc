// Copyright 2023 Christoph Weis
#include "include/materials.h"

#include <algorithm>
#include <functional>

#include "include/colors.h"
#include "include/math.h"
#include "materials.h"

Vec3 LambertianInverseRay(Vec3 object_normal) {
  return -(object_normal + RandomUnitVector()).NormalizedNonzero();
}

ScatterData Material::InverseScatter(
    const ReferenceFrameHit &rf_hit,
    SpectrumTransform &cumulative_transform) const {
  if (dielectric_) {
    return DielectricInverseScatter(rf_hit, cumulative_transform);
  }
  float scatter_selector = RandomReal() * non_dielectric_multipliers_();

  if (scatter_selector < diffuse_) {
    return LambertianInverseScatter(rf_hit, cumulative_transform);
  } else if (scatter_selector < diffuse_ + specular_) {
    return SpecularInverseScatter(rf_hit, cumulative_transform);
  } else {
    return SubsurfaceInverseScatter(rf_hit, cumulative_transform);
  }
}
Material &Material::SetProperty(std::string scalar_property_name,
                                float new_value) {
  if (scalar_property_name == "albedo") {
    albedo_ = new_value;
  } else if (scalar_property_name == "diffuse") {
    diffuse_ = new_value;
  } else if (scalar_property_name == "specular") {
    specular_ = new_value;
  } else if (scalar_property_name == "fuzz") {
    fuzz_ = new_value;
  } else if (scalar_property_name == "dielectric_n") {
    dielectric_n_ = new_value;
  } else if (scalar_property_name == "subsurface_scatters") {
    subsurface_scatters_ = new_value;
  } else {
    throw std::invalid_argument("There's no scalar property with name '" +
                                scalar_property_name + "'.");
  }
  return *this;
}

Material &Material::SetProperty(std::string spectral_property_name,
                                Spectrum new_value) {
  if (spectral_property_name == "emission") {
    emission_ = new_value;
  } else if (spectral_property_name == "absorption") {
    absorption_ = new_value;
  } else if (spectral_property_name == "subsurface absorption") {
    subsurface_absorption_ = new_value;
  } else {
    throw std::invalid_argument("There's no spectral property with name '" +
                                spectral_property_name + "'.");
  }
  return *this;
}

Material &Material::SetProperty(std::string bool_property_name,
                                bool new_value) {
  if (bool_property_name == "dielectric") {
    dielectric_ = new_value;
  } else {
    throw std::invalid_argument("There's no boolean property with name '" +
                                bool_property_name + "'.");
  }
  return *this;
}

Spectrum Material::EmissionSpectrum(const ReferenceFrameHit &rf_hit) const {
  float dot_factor = Dot3(rf_hit.normal, rf_hit.scattered);
  if (dot_factor < 0) {
    return kBlack;
  }  // else
  return Dot3(rf_hit.normal, rf_hit.scattered) * emission_;
}

Spectrum Material::AbsorptionCurve(const ReferenceFrameHit &rf_hit) const {
  // No dot product factor here, because the absorption curve specifies
  // the *fraction* of light rays absorbed depending on wavelength
  return absorption_;
}

ScatterData Material::LambertianInverseScatter(
    const ReferenceFrameHit &rf_hit,
    SpectrumTransform &cumulative_transform) const {
  cumulative_transform.ApplyFactor(albedo_ * non_dielectric_multipliers_());
  cumulative_transform.ApplyAbsorption(absorption_);
  return ScatterData(LambertianInverseRay(rf_hit.normal));
}

ScatterData Material::SpecularInverseScatter(
    const ReferenceFrameHit &rf_hit,
    SpectrumTransform &cumulative_transform) const {
  cumulative_transform.ApplyFactor(albedo_ * non_dielectric_multipliers_());
  cumulative_transform.ApplyAbsorption(absorption_);

  Vec3 reflected_vel = rf_hit.scattered.Reflect(rf_hit.normal);
  // apply fuzz
  return ScatterData(
      (reflected_vel.NormalizedNonzero() + fuzz_ * RandomVectorInUnitBall())
          .NormalizedNonzero() *
      reflected_vel.Norm());
}

ScatterData Material::DielectricInverseScatter(
    const ReferenceFrameHit &rf_hit,
    SpectrumTransform &cumulative_transform) const {
  cumulative_transform.ApplyFactor(albedo_);
  cumulative_transform.ApplyAbsorption(absorption_);
  // Adapted from “Ray Tracing in One Weekend.”
  // raytracing.github.io/books/RayTracingInOneWeekend.html
  // (we need to keep track of the group velocity of the ray.)
  Vec3 ray_vel = rf_hit.scattered;
  Vec3 normal = rf_hit.normal;

  float n_ratio;
  if (rf_hit.outside) {
    n_ratio = 1 / dielectric_n_;
  } else {
    n_ratio = dielectric_n_;
  }
  float n_ratio_2 = n_ratio * n_ratio;

  float cos_theta = ray_vel.NormalizedNonzero().Dot3(normal);
  float new_sin_theta = sqrtf(1 - cos_theta * cos_theta) * n_ratio;

  if (new_sin_theta > 1.0f) {
    // reflect fully
    return ScatterData(ray_vel.Reflect(normal));
  }

  // Schlick's approximation
  float r0 = (1 - n_ratio) / (1 + n_ratio);
  r0 = r0 * r0;
  float reflectance = r0 + (1 - r0) * powf(1 - cos_theta, 5);

  float reflect_or_refract = RandomReal();

  if (reflect_or_refract < reflectance) {
    return ray_vel.Reflect(normal);
  }  // else

  Vec3 parallel = Dot3(ray_vel, normal) * normal;
  Vec3 new_perpendicular = n_ratio_2 * (ray_vel - parallel);
  Vec3 new_parallel =
      sqrtf(ray_vel.NormSq() * n_ratio_2 - new_perpendicular.NormSq()) * normal;
  Vec3 refracted_vel = new_perpendicular + new_parallel;
  return ScatterData(refracted_vel);
}

ScatterData Material::SubsurfaceInverseScatter(
    const ReferenceFrameHit &rf_hit,
    SpectrumTransform &cumulative_transform) const {
  cumulative_transform.ApplyFactor(albedo_ * non_dielectric_multipliers_());
  cumulative_transform.ApplyAbsorption(subsurface_absorption_);

  return ScatterData(LambertianInverseRay(rf_hit.normal));
}

Material DielectricMaterial(float dielectric_n) {
  return Material(1.0f, 0.0f, 0.0f, 0.0f, kBlack, kBlack, true, dielectric_n);
}

Material Metal(float albedo, Spectrum absorption, float fuzz) {
  return Material(albedo, 0.0f, 1.0f, fuzz, absorption);
}

Material Mirror(float fuzz) { return Metal(1.0f, kBlack, fuzz); }

Material Light(Spectrum emission, Spectrum absorption) {
  return Metal(0.0f, absorption, 0.0f).SetProperty("emission", emission);
}
