// Copyright 2023 Christoph Weis
#include "include/colors.h"

#include <algorithm>
#include <cmath>
#include <functional>

#include "colors.h"
#include "include/math.h"

/*
  ABSORPTION THINGS
*/
const float kAbsorptionThreshold =
    1e-2;  // disregard absorption modes with area below this threshold

/*
  GAUSSIAN THINGS
*/

const float kStdGaussianPrefactor = 0.3989422804f;       // 1/sqrt(2 pi)
const float kInvStdGaussianPrefactor = 2.506628274631f;  // sqrt(2 pi)

// use non-positive standard deviation as code for infinite standard deviation
inline float piecewise_gaussian(float x, float a, float mu, float sigma_left,
                                float sigma_right) {
  if (x <= mu) {
    if (sigma_left <= 0) {
      return a;
    }
    return a * expf(-(x - mu) * (x - mu) / (2 * sigma_left * sigma_left));
  }
  if (sigma_right <= 0) {
    return a;
  }
  return a * expf(-(x - mu) * (x - mu) / (2 * sigma_right * sigma_right));
}

inline float std_gaussian(float x) {
  return kStdGaussianPrefactor * expf(-x * x / 2);
}

const float kGaussianIntegral_p = 0.33267f;
const float kGaussianIntegral_a1 = 0.4361836f;
const float kGaussianIntegral_a2 = -0.1201676f;
const float kGaussianIntegral_a3 = 0.9372980f;

inline float std_gaussian_integral(float x) {
  // After Zelen & Severo (1964) algorithm 26.2.16 (error ~ 1e-5)
  if (x < 0) {
    return 0.5f - std_gaussian_integral(-x);
  }
  float t = 1 / (1 + kGaussianIntegral_p * x);
  return 1 - std_gaussian(x) * t *
                 (kGaussianIntegral_a1 +
                  t * (kGaussianIntegral_a2 + t * kGaussianIntegral_a3));
  /* (fewer multiplications than the below)
  return 1 - std_gaussian(x) *
                 ((kGaussianIntegral_a1 * t) + (kGaussianIntegral_a2 * t * t) +
                  (kGaussianIntegral_a3 * t * t * t));
  */
}

/*
const float kGaussianIntegral_p = 0.2316419f;
const float kGaussianIntegral_b1 = 0.319381530f;
const float kGaussianIntegral_b2 = -0.356563782f;
const float kGaussianIntegral_b3 = 1.781477937f;
const float kGaussianIntegral_b4 = -1.821255978f;
const float kGaussianIntegral_b5 = 1.330274429f;

inline
float std_gaussian_integral(float x) {
// After Zelen & Severo (1964) algorithm 26.2.17 (error ~ 1e-7)
  if (x < 0) {
    return 0.5f - std_gaussian_integral(-x);
  }
  float t = 1/(1+kGaussianIntegral_p * x);
  // return 1 - std_gaussian(x) * (
  // (kGaussianIntegral_b1 * t) +
  // (kGaussianIntegral_b2 * t * t) +
  // (kGaussianIntegral_b3 * t * t * t) +
  // (kGaussianIntegral_b4 * t * t * t * t) +
  // (kGaussianIntegral_b5 * t * t * t * t * t)
  // );
  return 1 - std_gaussian(x) * t * (
  kGaussianIntegral_b1 + t * (
  kGaussianIntegral_b2 + t * (
  kGaussianIntegral_b3 + t * (
  kGaussianIntegral_b4 + t *
  kGaussianIntegral_b5))));
}
*/

float gaussian_integral(float x, float a, float mu, float sigma) {
  return a * sigma * std_gaussian_integral((x - mu) / sigma);
}

inline float product_mu(float mu1, float sig1, float mu2, float sig2) {
  if (sig1 >= 0 && sig2 >= 0) {
    float var1 = sig1 * sig1;
    float var2 = sig2 * sig2;
    return (var2 * mu1 + var1 * mu2) / (var1 + var2);
  } else if (sig1 >= 0) {
    return mu1;
  } else {
    return mu2;
  }
}

inline float product_sig(float sig1, float sig2) {
  if (sig1 >= 0 && sig2 >= 0) {
    float var1 = sig1 * sig1;
    float var2 = sig2 * sig2;
    return sqrtf((var1 * var2) / (var1 + var2));
  } else if (sig1 >= 0) {
    return sig1;
  } else {
    return sig2;
  }
}

inline float product_amp_factor(float mu1, float sig1, float mu2, float sig2) {
  if (sig1 >= 0 && sig2 >= 0) {
    return expf(-(mu1 - mu2) * (mu1 - mu2) / (2 * (sig1 * sig1 + sig2 * sig2)));
  }
  return 1;
}

/*
  Color formats
*/

rgbData operator+(const rgbData &left, const rgbData &right) {
  return rgbData(left.r + right.r, left.g + right.g, left.b + right.b);
}

rgbData operator*(const rgbData &rgb, float factor) {
  return rgbData(rgb.r * factor, rgb.g * factor, rgb.b * factor);
}

rgbData operator*(float factor, const rgbData &rgb) { return rgb * factor; }

rgbData operator/(const rgbData &rgb, float div) { return rgb * (1 / div); }

XYZData operator+(const XYZData &left, const XYZData &right) {
  return XYZData(left.X + right.X, left.Y + right.Y, left.Z + right.Z);
}

XYZData operator-(const XYZData &left, const XYZData &right) {
  return XYZData(left.X - right.X, left.Y - right.Y, left.Z - right.Z);
}

XYZData operator*(const XYZData &XYZ, float factor) {
  return XYZData(XYZ.X * factor, XYZ.Y * factor, XYZ.Z * factor);
}

XYZData operator*(float factor, const XYZData &XYZ) { return XYZ * factor; }

XYZData operator/(const XYZData &XYZ, float div) {
  return XYZData(XYZ.X / div, XYZ.Y / div, XYZ.Z / div);
}

// ----------------------
//  COLOR CONVERSIONS
// ----------------------

/*
  XYZ -> RGB
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
  XYZ -> rgb
*/
XYZData &XYZData::operator+=(const XYZData &other) {
  X += other.X;
  Y += other.Y;
  Z += other.Z;
  return *this;
}
XYZData &XYZData::operator-=(const XYZData &other) {
  X -= other.X;
  Y -= other.Y;
  Z -= other.Z;
  return *this;
}
XYZData &XYZData::operator*=(float factor) {
  X *= factor;
  Y *= factor;
  Z *= factor;
  return *this;
}
XYZData &XYZData::operator/=(float div) {
  X /= div;
  Y /= div;
  Z /= div;
  return *this;
}

rgbData XYZData::To_rgb() const {
  return rgbData(2.36461385f * X - 0.89654057f * Y - 0.46807328f * Z,
                 -0.51516621f * X + 1.4264081f * Y + 0.0887581f * Z,
                 0.0052037f * X - 0.01440816f * Y + 1.00920446f * Z);
}

/*
  rgb -> RGB
*/
RGBData rgbData::ToRGB(float rescale_factor, int max) const {
  return RGBData{.R = std::clamp(static_cast<int>(rescale_factor * r), 0, max),
                 .G = std::clamp(static_cast<int>(rescale_factor * g), 0, max),
                 .B = std::clamp(static_cast<int>(rescale_factor * b), 0, max)};
}

rgbData &rgbData::operator+=(const rgbData &other) {
  r += other.r;
  g += other.g;
  b += other.b;
  return *this;
}

rgbData &rgbData::operator-=(const rgbData &other) {
  r -= other.r;
  g -= other.g;
  b -= other.b;
  return *this;
}

rgbData &rgbData::operator*=(float factor) {
  r *= factor;
  g *= factor;
  b *= factor;
  return *this;
}

rgbData &rgbData::operator/=(float div) {
  r /= div;
  g /= div;
  b /= div;
  return *this;
}

std::vector<GaussianData> Spectrum::getModes() const {
  std::vector<GaussianData> modes(modes_);

  for (GaussianData &gaussian : modes) {
    gaussian.Rescale(amplitude_factor_, scale_factor_);
  }

  return modes;
}

size_t Spectrum::size() const { return modes_.size(); }

Spectrum &Spectrum::operator+=(const Spectrum &other) {
  for (const GaussianData &mode : other.getModes()) {
    modes_.push_back(mode.Rescaled(other.amplitude_factor_ / amplitude_factor_,
                                   other.scale_factor_ / scale_factor_));
  }
  return *this;
}

Spectrum Spectrum::operator+(const Spectrum &other) const {
  std::vector<GaussianData> modes(getModes());
  std::vector<GaussianData> other_modes(other.getModes());
  modes.insert(modes.end(), other_modes.cbegin(), other_modes.cend());
  return Spectrum(modes);
}

Spectrum &Spectrum::operator+=(const GaussianData &other) {
  modes_.push_back(
      other.Rescaled(1.0f / amplitude_factor_, 1.0f / scale_factor_));
  return *this;
}

Spectrum Spectrum::operator+(const GaussianData &other) const {
  Spectrum sum(*this);
  sum += other;
  return sum;
}

Spectrum &Spectrum::operator*=(float factor) {
  amplitude_factor_ *= factor;
  return *this;
}

Spectrum Spectrum::operator*(float factor) const {
  Spectrum product(*this);
  product *= factor;
  return product;
}

Spectrum operator*(float factor, const Spectrum &spec) { return spec * factor; }

/*
  Black body spectrum fit with three Gaussian Modes
  (see python script 'fitBBrad.py' in 'misc' subfolder).
*/
GaussianData kBlackBodyMode1 = GaussianData(2.582f, 0.2832f, 0.08492f);
GaussianData kBlackBodyMode2 = GaussianData(1.928f, 0.4818f, 0.1804f);
GaussianData kBlackBodyMode3 = GaussianData(1.751f, 0.1812f, 0.04919f);

float kBlackBodyRescale = 14387768.8;  // h*c/K in units [nm K]

// TODO(c): figure out normalising factor for area
Spectrum BlackBodyRadiation(float T, float amplitude_rescale) {
  // rescale appropriately
  return Spectrum({kBlackBodyMode1, kBlackBodyMode2, kBlackBodyMode3})
      .Rescaled(kBlackBodyRescale / T, amplitude_rescale);
}

Spectrum &Spectrum::Rescale(float amplitude_factor, float scale_factor) {
  for (GaussianData &mode : modes_) {
    mode.Rescale(amplitude_factor, scale_factor);
  }
  return *this;
}

Spectrum Spectrum::Rescaled(float amplitude_factor, float scale_factor) const {
  std::vector<GaussianData> modes(getModes());
  for (GaussianData &mode : modes) {
    mode.Rescale(amplitude_factor, scale_factor);
  }
  return Spectrum(modes);
}

SpectrumTransform &SpectrumTransform::ApplyFactor(float factor) {
  brightness_ *= factor;
  return *this;
}

rgbData SpectrumTransform::ColorFrom(const Spectrum &spectrum) const {
  // Ultimately, we want to evaluate
  // emission * (1-A1) * (1-A2) * ... =
  // emission * (1-A)
  // where A = (A1+A2+...) - (A1*A2+A1*A3+...) + ...
  // We do this by recursively evaluating this product,
  // adding in A1, then A2, then A3, ...
  // The products calculated along the way are stored
  // in the iterative_absorption_mode_buffer.

  // Create absorption buffer representing A
  XYZData color;
  std::vector<GaussianData> iterative_absorption_mode_buffer;
  for (const Spectrum &absorption : camera_frame_absorptions_) {
    std::vector<GaussianData> old_modes(iterative_absorption_mode_buffer);

    for (const GaussianData &new_absorption_mode : absorption.getModes()) {
      // speed-up: skip low-absorption modes
      if (abs(new_absorption_mode.amplitude * new_absorption_mode.sig) <
          kAbsorptionThreshold)
        continue;

      iterative_absorption_mode_buffer.push_back(new_absorption_mode);
      for (const GaussianData &old_mode : old_modes) {
        GaussianData mode =
            old_mode * new_absorption_mode.Rescaled(-1.0f, 1.0f);
        if (abs(mode.amplitude * mode.sig < kAbsorptionThreshold)) continue;
        iterative_absorption_mode_buffer.push_back(mode);
      }
    }
  }

  // Rescale emission according to transform, then compute resulting color.
  for (const GaussianData &em_mode :
       spectrum.Rescaled(brightness_, wavelength_rescale_).getModes()) {
    color += em_mode.Color();
    for (const GaussianData &ab_mode : iterative_absorption_mode_buffer) {
      color -= (em_mode * ab_mode).Color();
    }
  }

  return color.To_rgb();
}

SpectrumTransform &SpectrumTransform::ApplyAbsorption(const Spectrum &absorb) {
  // Add the absorption spectrum, transformed to camera frame.
  // Note it specifies a *fraction* as a function of wavelength,
  // so its amplitude transforms trivially.
  if (absorb.size() == 0) {
    return *this;
  }  // else
  camera_frame_absorptions_.push_back(
      absorb.Rescaled(1.0f, wavelength_rescale_));

  return *this;
}

SpectrumTransform &SpectrumTransform::ApplyTransformationToFrame(
    const Vec3 &ray_vel, const Vec3 &rf_vel) {
  // Both velocities must be specified in standard frame
  // Apply the Doppler rescaling of frequencies and spectral radiance
  float doppler =
      rf_vel.Gamma() * (1 - Dot3(ray_vel.NormalizedNonzero(), rf_vel));
  float doppler_5 = doppler * doppler * doppler * doppler * doppler;
  // Specific radiative intensity (intensity per unit wavelength)
  // rescales by doppler**5, while wavelength scales by 1/doppler
  brightness_ *= doppler_5;
  wavelength_rescale_ /= doppler;
  return *this;
}

SpectrumTransform &SpectrumTransform::ApplyTransformationFromFrame(
    const Vec3 &ray_vel, const Vec3 &rf_vel) {
  // Both velocities must be specified in standard frame
  // Apply the Doppler rescaling of frequencies and spectral radiance
  // This is inverse to TransformationToFrame
  float doppler =
      rf_vel.Gamma() * (1 - Dot3(ray_vel.NormalizedNonzero(), rf_vel));
  float doppler_5 = doppler * doppler * doppler * doppler * doppler;
  // Specific radiative intensity (intensity per unit wavelength)
  // rescales by doppler**5, while wavelength scales by 1/doppler
  brightness_ /= doppler_5;
  wavelength_rescale_ *= doppler;
  return *this;
}

GaussianData GaussianData::Rescaled(float amplitude_factor,
                                    float scale_factor) const {
  return GaussianData(amplitude * amplitude_factor, mean * scale_factor,
                      sig * scale_factor);
}

GaussianData &GaussianData::Rescale(float amplitude_factor,
                                    float scale_factor) {
  amplitude *= amplitude_factor;
  mean *= scale_factor;
  sig *= scale_factor;
  return *this;
}

float GaussianData::Area() {
  return amplitude * sig * kInvStdGaussianPrefactor;
}

GaussianData GaussianData::operator*(float factor) const {
  float new_a = amplitude * factor;
  return GaussianData(new_a, mean, sig);
}

GaussianData &GaussianData::operator*=(float factor) {
  amplitude *= factor;
  return *this;
}

// The product of two Gaussians is a Gaussian
GaussianData GaussianData::operator*(const GaussianData &other) const {
  float new_a = amplitude * other.amplitude *
                product_amp_factor(mean, sig, other.mean, other.sig);
  float new_mean = product_mu(mean, sig, other.mean, other.sig);
  float new_sig = product_sig(sig, other.sig);

  return GaussianData(new_a, new_mean, new_sig);
}

GaussianData &GaussianData::operator*=(const GaussianData &other) {
  amplitude *=
      other.amplitude * product_amp_factor(mean, sig, other.mean, other.sig);
  mean = product_mu(mean, sig, other.mean, other.sig);
  sig = product_sig(sig, other.sig);
  return *this;
}

// TODO(c): namespacing
/*
    Compute XYZData from wavelength following
    Wyman, Chris, Peter-Pike Sloan, and Peter Shirley.
    "Simple analytic approximations to the CIE XYZ color
    matching functions." J. Comput. Graph. Tech 2.2 (2013): 11.
  */
const GaussianData kCIEx1 = GaussianData(1.042f, 596.4f, 34.34f);
const GaussianData kCIEx2 = GaussianData(0.3687f, 446.7f, 19.50f);
const GaussianData kCIEy = GaussianData(1.022f, 559.1f, 41.80f);
const GaussianData kCIEz = GaussianData(1.859f, 451.0f, 22.74f);

XYZData GaussianData::Color() const {
  const GaussianData &mode = *this;
  return XYZData((kCIEx1 * mode).Area() + (kCIEx2 * mode).Area(),
                 (kCIEy * mode).Area(), (kCIEz * mode).Area());
}
