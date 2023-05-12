// Copyright 2023 Christoph Weis
#ifndef COLORS_H_INCLUDED
#define COLORS_H_INCLUDED

#include <vector>

#include "math.h"

const float kPi = 3.141592653589793f;

/*
  COLOR STRUCTURES
*/
struct RGBData {
  int R, G, B;
};

struct rgbData {
  float r, g, b;
  rgbData(float r_ = 0, float g_ = 0, float b_ = 0) : r(r_), g(g_), b(b_) {}
  RGBData ToRGB(float rescale_factor = 255.0f, int max = 255) const;
  rgbData &operator+=(const rgbData &);
  rgbData &operator-=(const rgbData &);
  rgbData &operator*=(float);
  rgbData &operator/=(float);
};
rgbData operator+(const rgbData &, const rgbData &);
rgbData operator*(const rgbData &, float);
rgbData operator*(float, const rgbData &);
rgbData operator/(const rgbData &, float);

struct XYZData {
  float X, Y, Z;
  // compute XYZ from wavelength
  XYZData(float _X = 0, float _Y = 0, float _Z = 0) : X(_X), Y(_Y), Z(_Z) {}

  XYZData &operator+=(const XYZData &);
  XYZData &operator-=(const XYZData &);
  XYZData &operator*=(float);
  XYZData &operator/=(float);

  rgbData To_rgb() const;
};
XYZData operator+(const XYZData &, const XYZData &);
XYZData operator-(const XYZData &, const XYZData &);
XYZData operator*(const XYZData &, float);
XYZData operator*(float, const XYZData &);
XYZData operator/(const XYZData &, float);

// GaussianData contains parameters of a Gaussian function
struct GaussianData {
  GaussianData(float _amplitude, float _mean, float _sig)
      : amplitude(_amplitude), mean(_mean), sig(_sig) {}

  GaussianData operator*(float factor) const;
  GaussianData operator*(const GaussianData &other) const;
  GaussianData &operator*=(float factor);
  GaussianData &operator*=(const GaussianData &other);
  GaussianData Rescaled(float amplitude_factor, float scale_factor) const;
  GaussianData &Rescale(float amplitude_factor, float scale_factor);
  float Area();
  XYZData Color() const;

  float amplitude;
  float mean;
  float sig;
};

class Spectrum {
 public:
  explicit Spectrum(std::vector<GaussianData> modes = {})
      : amplitude_factor_(1.0f), scale_factor_(1.0f), modes_(modes) {}
  explicit Spectrum(GaussianData mode)
      : Spectrum(std::vector<GaussianData>{mode}) {}
  std::vector<GaussianData> getModes() const;

  Spectrum &operator+=(const Spectrum &other);
  Spectrum operator+(const Spectrum &other) const;
  Spectrum &operator+=(const GaussianData &other);
  Spectrum operator+(const GaussianData &other) const;
  Spectrum &operator*=(float factor);
  Spectrum operator*(float factor) const;

  Spectrum Rescaled(float amplitude_factor, float scale_factor) const;
  Spectrum &Rescale(float amplitude_factor, float scale_factor);

 private:
  float amplitude_factor_;  // multiplies amplitudes
  float scale_factor_;      // rescales wavelengths
  std::vector<GaussianData> modes_;
};
Spectrum operator*(float factor, const Spectrum &spec);

// BLACKBODY SPECTRUM
Spectrum BlackBodyRadiation(float T, float amplitude_rescale = 1.0f);

class SpectrumTransform {
 public:
  SpectrumTransform(float brightness = 1.0f, float wavelength_rescale = 1.0f,
                    std::vector<Spectrum> absorbtions = {})
      : brightness_(brightness),
        wavelength_rescale_(wavelength_rescale),
        camera_frame_absorbtions_(absorbtions) {}
  SpectrumTransform(float brightness, float wavelength_rescale, Spectrum absorb)
      : SpectrumTransform(brightness, wavelength_rescale,
                          std::vector<Spectrum>{absorb}) {}

  rgbData ColorFrom(const Spectrum &spectrum) const;
  SpectrumTransform &ApplyAbsorption(const Spectrum &absorb_);
  SpectrumTransform &ApplyFactor(float factor);
  SpectrumTransform &ApplyTransformationToFrame(const Vec3 &ray_vel,
                                                const Vec3 &rf_vel);
  SpectrumTransform &ApplyTransformationFromFrame(const Vec3 &ray_vel,
                                                  const Vec3 &rf_vel);

 private:
  float brightness_;  // overall multiplier

  float wavelength_rescale_;  // rescale of wavelengths in spectrum
  std::vector<Spectrum>
      camera_frame_absorbtions_;  // absorption spectrum to be applied
};

const Spectrum kBlack;
const Spectrum kWhite(GaussianData(1.0f, 500.0f, -1.0f));

const float kNarrowBandSig = 0.5f;
const float kNarrowBandBrightAmplitude = 1.0f;
const float kNarrowBandLightAmplitude = 0.5f;
const float kNarrowBandDimAmplitude = 0.1f;

inline Spectrum BrightNarrowBand(float lambda) {
  return Spectrum(
      GaussianData(kNarrowBandBrightAmplitude, lambda, kNarrowBandSig));
}

inline Spectrum LightNarrowBand(float lambda) {
  return Spectrum(
      GaussianData(kNarrowBandBrightAmplitude, lambda, kNarrowBandSig));
}

inline Spectrum DimNarrowBand(float lambda) {
  return Spectrum(
      GaussianData(kNarrowBandBrightAmplitude, lambda, kNarrowBandSig));
}

#endif