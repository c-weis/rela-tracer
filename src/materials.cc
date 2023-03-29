// Copyright 2023 Christoph Weis
#include "include/materials.h"

#include <algorithm>
#include <functional>

#include "include/math.h"

inline float piecewise_gaussian(float x, float mu, float sigma_left,
                                float sigma_right) {
  if (x <= mu) {
    return expf(-(x - mu) * (x - mu) / (2 * sigma_left * sigma_left));
  }
  return expf(-(x - mu) * (x - mu) / (2 * sigma_right * sigma_right));
}

rgbData operator*(const rgbData &rgb, float factor) {
  return rgbData{
      .r = rgb.r * factor,
      .g = rgb.g * factor,
      .b = rgb.b * factor,
  };
}

rgbData operator*(float factor, const rgbData &rgb) { return rgb * factor; }

rgbData operator/(const rgbData &rgb, float div) { return rgb * (1 / div); }

XYZData operator+(const XYZData &left, const XYZData &right) {
  return XYZData(left.X + right.X, left.Y + right.Y, left.Z + right.Z);
}

XYZData operator*(const XYZData &XYZ, float factor) {
  return XYZData(XYZ.X * factor, XYZ.Y * factor, XYZ.Z * factor);
}

XYZData operator*(float factor, const XYZData &XYZ) { return XYZ * factor; }

XYZData operator/(const XYZData &XYZ, float div) {
  return XYZData(XYZ.X / div, XYZ.Y / div, XYZ.Z / div);
}

/*
    Compute XYZData from wavelength following
    Wyman, Chris, Peter-Pike Sloan, and Peter Shirley.
    "Simple analytic approximations to the CIE XYZ color
    matching functions." J. Comput. Graph. Tech 2.2 (2013): 11.
  */
XYZData::XYZData(float lambda)
    : X(1.056 * piecewise_gaussian(lambda, 599.8, 37.9, 31.0) +
        0.362 * piecewise_gaussian(lambda, 442.0, 16.0, 26.7) -
        0.065 * piecewise_gaussian(lambda, 501.1, 20.4, 26.2)),
      Y(0.821 * piecewise_gaussian(lambda, 568.8, 46.9, 40.5) +
        0.286 * piecewise_gaussian(lambda, 530.9, 16.3, 31.1)),
      Z(1.217 * piecewise_gaussian(lambda, 437.0, 11.8, 36.0) +
        0.681 * piecewise_gaussian(lambda, 459.0, 26.0, 13.8)){};

/*
  Convert XYZ to RGB using standard linear map.
*/
RGBData RGBFromXYZ(XYZData XYZ) {
  auto [X, Y, Z] = XYZ;

  return {.R = static_cast<int>(
              255 * (2.36461385f * X - 0.89654057f * Y - 0.46807328f * Z)),
          .G = static_cast<int>(
              255 * (-0.51516621f * X + 1.4264081f * Y + 0.0887581f * Z)),
          .B = static_cast<int>(
              255 * (0.0052037f * X - 0.01440816f * Y + 1.00920446f * Z))};
}

/*
  Convert XYZ to rgb using standard linear map.
*/
rgbData XYZData::To_rgb() const {
  return {.r = 2.36461385f * X - 0.89654057f * Y - 0.46807328f * Z,
          .g = -0.51516621f * X + 1.4264081f * Y + 0.0887581f * Z,
          .b = 0.0052037f * X - 0.01440816f * Y + 1.00920446f * Z};
}

/*
  Convert ColorData to rgb data
*/
rgbData ColorData::To_rgb() const {
  XYZData XYZ(0, 0, 0);
  for (LightRay light_ray : light_rays) {
    float lambda =
        2 * M_PI / light_ray.k.t;  // the 0th ("time") coordinate of a
                                   // wave-4-vector is the angular velocity
                                   // omega (recall we set c=1)
    XYZ = XYZ + XYZData(lambda) * light_ray.brightness;
  }
  return XYZ.To_rgb();
}

/*
  Convert ColorData to RGB data
*/
RGBData ColorData::ToRGB() const {
  XYZData XYZ(0, 0, 0);
  for (LightRay light_ray : light_rays) {
    float lambda =
        2 * M_PI / light_ray.k.t;  // the 0th ("time") coordinate of a
                                   // wave-4-vector is the angular velocity
                                   // omega (recall we set c=1)
    XYZ = XYZ + XYZData(lambda) * light_ray.brightness;
  }
  return RGBFromXYZ(XYZ);
}

/*
  The method below is very inefficient and unnecessary for our application.
  It is commented out to avoid accidental use. We only use +=.

ColorData operator+(const ColorData &left, const ColorData &right) {
  std::vector<LightRay> combined_rays;
  combined_rays.assign(left.light_rays.cbegin(), left.light_rays.cend());
  combined_rays.insert(combined_rays.end(), right.light_rays.cbegin(),
                       right.light_rays.cend());
  return ColorData(combined_rays);
}
*/

ColorData operator*(const ColorData &color, float factor) {
  std::vector<LightRay> new_rays;
  std::transform(color.light_rays.cbegin(), color.light_rays.cend(),
                 std::back_inserter(new_rays),
                 [factor](LightRay lr) { return lr * factor; });
  return ColorData(new_rays);
}

ColorData operator*(float factor, const ColorData &color) {
  return color * factor;
}

ColorData operator/(const ColorData &color, float divisor) {
  return color * (1 / divisor);
}

ColorData &ColorData::operator+=(const ColorData &other) {
  light_rays.insert(light_rays.end(), other.light_rays.cbegin(),
                    other.light_rays.cend());
  return *this;
}

ColorData &ColorData::operator*=(float factor) {
  std::for_each(light_rays.begin(), light_rays.end(),
                [factor](LightRay &light_ray) { light_ray *= factor; });
  return *this;
}

ColorData &ColorData::operator/=(float div) {
  std::for_each(light_rays.begin(), light_rays.end(),
                [div](LightRay &light_ray) { light_ray /= div; });
  return *this;
}

ColorData &ColorData::TransformToFrame(const Vec3 &frame_velocity) {
  std::for_each(light_rays.begin(), light_rays.end(),
                [frame_velocity](LightRay &light_ray) {
                  light_ray.TransformToFrame(frame_velocity);
                });
  return *this;
}

ColorData &ColorData::TransformFromFrame(const Vec3 &frame_velocity) {
  std::for_each(light_rays.begin(), light_rays.end(),
                [frame_velocity](LightRay &light_ray) {
                  light_ray.TransformFromFrame(frame_velocity);
                });
  return *this;
}

RGBData RGBData::gray(int shade) {
  return RGBData{.R = shade, .G = shade, .B = shade};
}

const ColorData MonochromaticLight::EmittedColor(
    const ReferenceFrameHit &rf_hit) const {
  Vec4 k(omega_, omega_ * rf_hit.scattered / rf_hit.scattered.NormSq());
  return ColorData(
      {LightRay(brightness_ * Dot3(rf_hit.scattered, rf_hit.normal), k)});
}

Vec3 LambertianScatter::RandomUnitVector(std::mt19937 *r_gen) {
  Vec3 vec(1, 1, 1);
  std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
  while (vec.NormSq() > 1 || vec.NormSq() == 0.0f) {
    vec = Vec3(dist(*r_gen), dist(*r_gen), dist(*r_gen));
  }
  return vec.NormalizedNonzero();
}

Vec3 LambertianScatter::LambertianInverseRay(Vec3 object_normal,
                                             std::mt19937 *r_gen) {
  return -(object_normal + LambertianScatter::RandomUnitVector(r_gen))
              .NormalizedNonzero();
}

std::vector<Vec3> LambertianScatter::InverseScatter(Vec3 scattered_ray,
                                                    Vec3 object_normal,
                                                    int max_number_of_rays,
                                                    std::mt19937 *r_gen) const {
  std::vector<Vec3> de_scatter_velocities;
  std::generate_n(std::back_inserter(de_scatter_velocities), max_number_of_rays,
                  [object_normal, r_gen]() {
                    return LambertianInverseRay(object_normal, r_gen);
                  });

  /* PASSED: Scatters anti-aligned with normal test
  std::for_each(de_scatter_velocities.cbegin(), de_scatter_velocities.cend(),
  [object_normal](Vec3 vel){ if (Dot3(object_normal, vel) >= 0) { std::cout <<
  "Error in Lambertian InverseScatter!" << std::endl;
    }
  });*/

  return de_scatter_velocities;
}

const ColorData LambertianScatter::ScatteredColor(
    const ReferenceFrameHit &rf_hit, const Vec3 &pre_scattered,
    const ColorData &pre_scatter_color) const {
  // TODO(c): implement absorption
  return pre_scatter_color;
}

std::vector<Vec3> Mirror::InverseScatter(
    Vec3 scattered_ray, Vec3 object_normal, int max_number_of_rays,
    std::mt19937 *random_number_generator) const {
  // reflect perfectly
  // TODO(c): change scatter code to allow weighting of rays
  //          then we can return a single ray here instead
  std::vector<Vec3> de_scatters;
  de_scatters.assign(max_number_of_rays,
                     scattered_ray.Reflected(object_normal));
  return de_scatters;
}

const ColorData Mirror::ScatteredColor(
    const ReferenceFrameHit &rf_hit, const Vec3 &pre_scattered,
    const ColorData &pre_scatter_color) const {
  // Mirror reflects color unchanged
  return pre_scatter_color;
}

LightRay operator*(const LightRay &light_ray, float factor) {
  return LightRay(light_ray.brightness * factor, light_ray.k);
}

LightRay operator*(float factor, const LightRay &light_ray) {
  return light_ray * factor;
}

LightRay &LightRay::operator*=(float factor) {
  brightness *= factor;
  return *this;
}

LightRay &LightRay::operator/=(float div) {
  brightness /= div;
  return *this;
}

LightRay LightRay::TransformedToFrame(const Vec3 &frame_velocity) const {
  float one_over_gamma = frame_velocity.GammaInv();
  return LightRay(brightness * one_over_gamma,
                  k.TransformedToFrame(frame_velocity));
}

LightRay LightRay::TransformedFromFrame(const Vec3 &frame_velocity) const {
  float gamma = frame_velocity.Gamma();
  return LightRay(brightness * gamma, k.TransformedFromFrame(frame_velocity));
}

LightRay &LightRay::TransformToFrame(const Vec3 &frame_velocity) {
  float one_over_gamma = frame_velocity.GammaInv();
  brightness *= one_over_gamma;
  k = k.TransformedToFrame(frame_velocity);
  return *this;
}

LightRay &LightRay::TransformFromFrame(const Vec3 &frame_velocity) {
  float gamma = frame_velocity.Gamma();
  brightness *= gamma;
  k = k.TransformedFromFrame(frame_velocity);
  return *this;
}

Material &Material::SetScatterer(Scatterer *scatterer) {
  scatterer_ = scatterer;
  return *this;
}

Material &Material::SetEmitter(Emitter *emitter) {
  emitter_ = emitter;
  return *this;
}

std::vector<Vec3> Material::InverseScatter(
    Vec3 scattered_ray, Vec3 object_normal, int max_number_of_rays,
    std::mt19937 *random_number_generator) const {
  if (scatterer_ == nullptr) {
    return std::vector<Vec3>();
  }
  return scatterer_->InverseScatter(scattered_ray, object_normal,
                                    max_number_of_rays,
                                    random_number_generator);
}

const ColorData Material::ScatteredColor(
    const ReferenceFrameHit &rf_hit, const Vec3 &pre_scattered,
    const ColorData &pre_scatter_color) const {
  if (scatterer_ == nullptr) {
    return kBlack;
  }
  return scatterer_->ScatteredColor(rf_hit, pre_scattered, pre_scatter_color);
}

const ColorData Material::EmittedColor(const ReferenceFrameHit &rf_hit) const {
  if (emitter_ == nullptr) {
    return kBlack;
  }
  return emitter_->EmittedColor(rf_hit);
}

RGBData rgbData::ToRGB(float rescale_factor, int max) const {
  return RGBData{.R = std::clamp(static_cast<int>(rescale_factor * r), 0, max),
                 .G = std::clamp(static_cast<int>(rescale_factor * g), 0, max),
                 .B = std::clamp(static_cast<int>(rescale_factor * b), 0, max)};
}
