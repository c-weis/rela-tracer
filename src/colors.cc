// Copyright 2023 Christoph Weis
#include "include/colors.h"

#include <algorithm>
#include <cmath>
#include <functional>

#include "colors.h"
#include "include/math.h"

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

Spectrum &Spectrum::operator+=(const Spectrum &other) {
  for (GaussianData mode : other.getModes()) {
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
  for (auto mode : modes_) {
    mode.Rescale(amplitude_factor, scale_factor);
  }
  return *this;
}

Spectrum Spectrum::Rescaled(float amplitude_factor, float scale_factor) const {
  std::vector<GaussianData> modes(getModes());
  for (auto mode : modes) {
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

  XYZData color;
  std::vector<GaussianData> iterative_absorption_mode_buffer;
  for (const Spectrum absorption : camera_frame_absorbtions_) {
    std::vector<GaussianData> old_modes(iterative_absorption_mode_buffer);

    for (const GaussianData new_absorption_mode : absorption.getModes()) {
      iterative_absorption_mode_buffer.push_back(new_absorption_mode);
      for (GaussianData old_mode : old_modes) {
        iterative_absorption_mode_buffer.push_back(
            old_mode * new_absorption_mode.Rescaled(-1.0f, 1.0f));
      }
    }
  }

  for (auto em_mode : spectrum.getModes()) {
    color += em_mode.Color();
    for (auto ab_mode : iterative_absorption_mode_buffer) {
      color -= (em_mode * ab_mode).Color();
    }
  }

  return color.To_rgb();
}

SpectrumTransform &SpectrumTransform::ApplyAbsorption(const Spectrum &absorb) {
  // Add the absorption spectrum, transformed to camera frame.
  // Note it specifies a *fraction* as a function of wavelength,
  // so its amplitude transforms trivially.
  camera_frame_absorbtions_.push_back(
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

// ---------------------------------------------
//   PIECEWISE GAUSSIAN THINGS
// ---------------------------------------------

/*
// PWGaussianData contains parameters necessary to describe
// a piecewise Gaussian function.
struct PWGaussianData {
  PWGaussianData(float _amplitude, float _mean, float _std_left,
                 float _std_right)
      : amplitude(_amplitude),
        mean(_mean),
        std_left(_std_left),
        std_right(_std_right) {}
  PWGaussianData(float amplitude, float mean, float std)
      : PWGaussianData(amplitude, mean, std, std) {}

  PWGaussianData Rescaled(float amplitude_factor, float scale_factor) const;
  PWGaussianData &Rescale(float amplitude_factor, float scale_factor);

  float amplitude;
  float mean;
  float std_left;
  float std_right;
};
XYZData EmittedColor(PWGaussianData emission);
XYZData AbsorbedColor(PWGaussianData emission, PWGaussianData absorption);
*/

/*
float piecewise_gaussian_integral(float x, float a, float mu, float sig_l,
                                  float sig_r) {
  if (x < mu) {
    return gaussian_integral(x, a, mu, sig_l);
  }
  return a * (sig_l - sig_r) * 0.5 + gaussian_integral(x, a, mu, sig_r);
}
*/

/* old version: compute via variance, doesn't implement negative=infinity hack
inline float product_mu(float mu1, float var1, float mu2, float var2) {
  return (var2 * mu1 + var1 * mu2) / (var1 + var2);
}

inline float product_var(float var1, float var2) {
  return (var1 * var2) / (var1 + var2);
}
*/

/*
    Compute XYZData from wavelength following
    Wyman, Chris, Peter-Pike Sloan, and Peter Shirley.
    "Simple analytic approximations to the CIE XYZ color
    matching functions." J. Comput. Graph. Tech 2.2 (2013): 11.
  */
// THESE ARE WRONG:
// const PWGaussianData kX1 = PWGaussianData(1.056f, 599.8f, 37.9f, 31.0f);
// const PWGaussianData kX2 = PWGaussianData(0.362f, 442.0f, 16.0f, 26.7f);
// const PWGaussianData kX3 = PWGaussianData(-0.065f, 501.1f, 20.4f, 26.2f);
// const PWGaussianData kY1 = PWGaussianData(0.821f, 568.8f, 46.9f, 40.5f);
// const PWGaussianData kY2 = PWGaussianData(0.286f, 530.9f, 16.3f, 31.1f);
// const PWGaussianData kZ1 = PWGaussianData(1.217f, 437.0f, 11.8f, 36.0f);
// const PWGaussianData kZ2 = PWGaussianData(0.681f, 459.0f, 26.0f, 13.8f);
//
// XYZData FromWavelength(float lambda) {
//   return XYZData((1.056 * piecewise_gaussian(lambda, 1.0f, 599.8, 37.9, 31.0)
//   +
//                   0.362 * piecewise_gaussian(lambda, 1.0f, 442.0, 16.0, 26.7)
//                   - 0.065 * piecewise_gaussian(lambda, 1.0f,
//                   501.1, 20.4, 26.2)),
//                  (0.821 * piecewise_gaussian(lambda, 1.0f, 568.8, 46.9, 40.5)
//                  +
//                   0.286 * piecewise_gaussian(lambda, 1.0f,
//                   530.9, 16.3, 31.1)),
//                  (1.217 * piecewise_gaussian(lambda, 1.0f, 437.0, 11.8, 36.0)
//                  +
//                   0.681 * piecewise_gaussian(lambda, 1.0f,
//                   459.0, 26.0, 13.8)));
// }
/*
float integrate_product(float a1, float mu1, float sig_l1, float sig_r1,
                        float a2, float mu2, float sig_l2, float sig_r2) {
  // make sure mu1 < mu2
  if (mu1 > mu2) {
    return integrate_product(a2, mu2, sig_l2, sig_r2, a1, mu1, sig_l1, sig_r1);
  }

  // The product of two piecewise Gaussians consists of three
  // Gaussians stitched together at mu1 and mu2.
  // Compute mu and sigma for the three sections,
  // then evaluate the integrals.
  float a = a1 * a2;

  float mu_left = product_mu(mu1, sig_l1, mu2, sig_l2);
  float sig_left = product_sig(sig_l1, sig_l2);
  float a_left = a * product_amp_factor(mu1, sig_l1, mu2, sig_l2);

  float mu_mid = product_mu(mu1, sig_r1, mu2, sig_l2);
  float sig_mid = product_sig(sig_r1, sig_l2);
  float a_mid = a * product_amp_factor(mu1, sig_r1, mu2, sig_l2);

  float mu_right = product_mu(mu1, sig_r1, mu2, sig_r2);
  float sig_right = product_sig(sig_r1, sig_r2);
  float a_right = a * product_amp_factor(mu1, sig_r1, mu2, sig_r2);

  float integral_left = gaussian_integral(mu1, a_left, mu_left, sig_left);
  float integral_mid = gaussian_integral(mu2, a_mid, mu_mid, sig_mid) -
                       gaussian_integral(mu1, a_mid, mu_mid, sig_mid);
  float integral_right = gaussian_integral(-mu2, a_right, mu_right, sig_right);

  return integral_left + integral_mid + integral_right;
}

float integrate_product(PWGaussianData gd1, PWGaussianData gd2) {
  return integrate_product(gd1.amplitude, gd1.mean, gd1.std_left, gd1.std_right,
                           gd2.amplitude, gd2.mean, gd2.std_left,
                           gd2.std_right);
}

float integrate_product(float a1, float mu1, float sig_l1, float sig_r1,
                        float a2, float mu2, float sig_l2, float sig_r2,
                        float a3, float mu3, float sig_l3, float sig_r3) {
  // make sure mu1 < mu2 < mu3 (bubble sort)
  if (mu2 > mu3) {
    return integrate_product(a1, mu1, sig_l1, sig_r1, a3, mu3, sig_l3, sig_r3,
                             a2, mu2, sig_l2, sig_r2);
  }
  if (mu1 > mu2) {
    return integrate_product(a2, mu2, sig_l2, sig_r2, a1, mu1, sig_l1, sig_r1,
                             a3, mu3, sig_l3, sig_r3);
  }

  // The product of three piecewise Gaussians consists of four
  // Gaussians stitched together at mu1, mu2 and mu3.
  // Compute mu and sigma for the four sections,
  // then evaluate the integrals.
  float a = a1 * a2 * a3;

  float mu12 = product_mu(mu1, sig_l1, mu2, sig_l2);
  float sig12 = product_sig(sig_l1, sig_l2);
  float amp12 = product_amp_factor(mu1, sig_l1, mu2, sig_l2);
  float mu_left = product_mu(mu12, sig12, mu3, sig_l3);
  float sig_left = product_sig(sig12, sig_l3);
  float a_left = a * amp12 * product_amp_factor(mu12, sig12, mu3, sig_l3);

  mu12 = product_mu(mu1, sig_r1, mu2, sig_l2);
  sig12 = product_sig(sig_r1, sig_l2);
  amp12 = product_amp_factor(mu1, sig_r1, mu2, sig_l2);
  float mu_midleft = product_mu(mu12, sig12, mu3, sig_l3);
  float sig_midleft = product_sig(sig12, sig_l3);
  float a_midleft = a * amp12 * product_amp_factor(mu12, sig12, mu3, sig_l3);

  mu12 = product_mu(mu1, sig_r1, mu2, sig_r2);
  sig12 = product_sig(sig_r1, sig_r2);
  amp12 = product_amp_factor(mu1, sig_r1, mu2, sig_r2);
  float mu_midright = product_mu(mu12, sig12, mu3, sig_l3);
  float sig_midright = product_sig(sig12, sig_l3);
  float a_midright = a * amp12 * product_amp_factor(mu12, sig12, mu3, sig_l3);

  float mu_right = product_mu(mu12, sig12, mu3, sig_r3);
  float sig_right = product_sig(sig12, sig_r3);
  float a_right = a * amp12 * product_amp_factor(mu12, sig12, mu3, sig_r3);

  float integral_left = gaussian_integral(mu1, a_left, mu_left, sig_left);
  float integral_midleft =
      gaussian_integral(mu2, a_midleft, mu_midleft, sig_midleft) -
      gaussian_integral(mu1, a_midleft, mu_midleft, sig_midleft);
  float integral_midright =
      gaussian_integral(mu3, a_midright, mu_midright, sig_midright) -
      gaussian_integral(mu2, a_midright, mu_midright, sig_midright);
  float integral_right = gaussian_integral(-mu3, a_right, mu_right, sig_right);

  return integral_left + integral_midleft + integral_midright + integral_right;
}

float integrate_product(PWGaussianData gd1, PWGaussianData gd2,
                        PWGaussianData gd3) {
  return integrate_product(gd1.amplitude, gd1.mean, gd1.std_left, gd1.std_right,
                           gd2.amplitude, gd2.mean, gd2.std_left, gd2.std_right,
                           gd3.amplitude, gd3.mean, gd3.std_left,
                           gd3.std_right);
}
*/

//
// XYZData EmittedColor(PWGaussianData emission) {
//   return XYZData(
//       integrate_product(kX1, emission) + integrate_product(kX2, emission) +
//           integrate_product(kX3, emission),
//       integrate_product(kY1, emission) + integrate_product(kY2, emission),
//       integrate_product(kZ1, emission) + integrate_product(kZ2, emission));
// }
//
// XYZData AbsorbedColor(PWGaussianData emission, PWGaussianData absorption) {
//   return XYZData(integrate_product(kX1, emission, absorption) +
//                      integrate_product(kX2, emission, absorption) +
//                      integrate_product(kX3, emission, absorption),
//                  integrate_product(kY1, emission, absorption) +
//                      integrate_product(kY2, emission, absorption),
//                  integrate_product(kZ1, emission, absorption) +
//                      integrate_product(kZ2, emission, absorption));
// }
