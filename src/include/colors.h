// Copyright 2023 Christoph Weis
#ifndef COLORS_H_INCLUDED
#define COLORS_H_INCLUDED

#include <vector>

#include "math.h"

const float kPi = 3.141592653589793f;

// -----------------------
//  COLOR DATA STRUCTURES
// -----------------------

// RGB values: R,G,B in [0,255]
struct RGBData {
  int R, G, B;
};

// rgb values: r,g,b in [0,1]
struct rgbData {
  float r, g, b;
  explicit rgbData(float r_ = 0, float g_ = 0, float b_ = 0)
      : r(r_), g(g_), b(b_) {}
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

// XYZ values: X,Y,Z reals
struct XYZData {
  float X, Y, Z;
  explicit XYZData(float _X = 0, float _Y = 0, float _Z = 0)
      : X(_X), Y(_Y), Z(_Z) {}

  XYZData &operator+=(const XYZData &);
  XYZData &operator-=(const XYZData &);
  XYZData &operator*=(float);
  XYZData &operator/=(float);

  // convert XYZ to rgb values
  rgbData To_rgb() const;
};
XYZData operator+(const XYZData &, const XYZData &);
XYZData operator-(const XYZData &, const XYZData &);
XYZData operator*(const XYZData &, float);
XYZData operator*(float, const XYZData &);
XYZData operator/(const XYZData &, float);

// ---------------
//  COLOR SPECTRA
// ---------------

// Parameters of a Gaussian function.
// Color spectra in our color model are sums of Gaussian modes.
struct Gaussian {
  // Gaussian mode with given amplitude, mean and std deviation (`_sig`).
  // sigma = -1 is code for a constant function
  Gaussian(float _amplitude, float _mean, float _sig)
      : amplitude(_amplitude), mean(_mean), sig(_sig) {}

  Gaussian operator*(float factor) const;
  Gaussian operator*(const Gaussian &other) const;
  Gaussian &operator*=(float factor);
  Gaussian &operator*=(const Gaussian &other);
  // Return rescaled Gaussian
  // Args:
  //  amplitude_factor: factor to apply to function value (usually brightness)
  //  scale_factor: factor to apply to coordinate (usually wavelength)
  Gaussian Rescaled(float amplitude_factor, float scale_factor) const;
  // Rescale this Gaussian
  // Args:
  //  amplitude_factor: factor to apply to function value (usually brightness)
  //  scale_factor: factor to apply to coordinate (usually wavelength)
  Gaussian &Rescale(float amplitude_factor, float scale_factor);
  // Return area of Gaussian.
  float Area();
  // Return XYZ color coordinates of this Gaussian considered as a mode
  // in a color spectrum. The x-coordinate is treated as the wavelength,
  // given in units of nanometers.
  XYZData Color() const;

  float amplitude;
  float mean;
  float sig;
};

// Color spectrum made up of Gaussian modes.
class Spectrum {
 public:
  // Spectrum with given Gaussian modes.
  explicit Spectrum(std::vector<Gaussian> modes = {})
      : amplitude_factor_(1.0f), scale_factor_(1.0f), modes_(modes) {}
  // Spectrum with single Gaussian mode.
  explicit Spectrum(Gaussian mode) : Spectrum(std::vector<Gaussian>{mode}) {}
  std::vector<Gaussian> getModes() const;
  // Returns number of Gaussian modes.
  size_t size() const;

  Spectrum &operator+=(const Spectrum &other);
  Spectrum operator+(const Spectrum &other) const;
  Spectrum &operator+=(const Gaussian &other);
  Spectrum operator+(const Gaussian &other) const;
  Spectrum &operator*=(float factor);
  Spectrum operator*(float factor) const;

  // Returns rescaled spectrum.
  // Args:
  //  amplitude_factor: factor rescaling brightness
  //  scale_factor: factor rescaling wavelength
  Spectrum Rescaled(float amplitude_factor, float scale_factor) const;
  // Rescales this spectrum.
  // Args:
  //  amplitude_factor: factor rescaling brightness
  //  scale_factor: factor rescaling wavelength
  Spectrum &Rescale(float amplitude_factor, float scale_factor);

  // Approximate Black Body Spectrum, rescaled such that the sun,
  // with temperature T=5800K has brightness amplitude 1.
  // Args:
  //  T: temperature in Kelvin
  //  amplitude_rescale: brightness rescale factor
  static Spectrum BlackBodyRadiation(float T, float amplitude_rescale = 1.0f);

  // Spectrum with amp=1, std deviation sig=0.5 nm, centred on given wavelength.
  static Spectrum NarrowBand(float lambda) {
    return Spectrum(Gaussian(1.0f, lambda, 0.5f));
  }

  // Empty spectrum
  static const Spectrum Black;
  // Fake white: constant across ALL wavelengths
  static const Spectrum White;

  // Narrowband Color representatives
  static const Spectrum Red, Orange, Yellow, Green, Cyan, Blue, Purple;

 private:
  float amplitude_factor_;       // multiplies amplitudes
  float scale_factor_;           // rescales wavelengths
  std::vector<Gaussian> modes_;  // contains spectral modes
};
Spectrum operator*(float factor, const Spectrum &spec);

// Transformations to be applied to color spectrum.
// Comprises brightness rescaling, frequency rescaling and absorptions.
class SpectrumTransform {
 public:
  SpectrumTransform(float brightness = 1.0f, float wavelength_rescale = 1.0f,
                    std::vector<Spectrum> absorptions = {})
      : brightness_(brightness),
        wavelength_rescale_(wavelength_rescale),
        camera_frame_absorptions_(absorptions) {}
  SpectrumTransform(float brightness, float wavelength_rescale, Spectrum absorb)
      : SpectrumTransform(brightness, wavelength_rescale,
                          std::vector<Spectrum>{absorb}) {}

  // Apply transform to given Spectrum and compute the resulting rgb color.
  rgbData ColorFrom(const Spectrum &spectrum) const;

  // Add absorption to transform.
  SpectrumTransform &ApplyAbsorption(const Spectrum &absorb_);
  // Add brightness factor to transform.
  SpectrumTransform &ApplyFactor(float factor);
  // Add Lorentz transformation into a frame to transform.
  // Args:
  //  ray_vel: velocity of light ray whose color spectrum is to be transformed
  //  rf_vel: velocity of reference frame into which we transform
  SpectrumTransform &ApplyTransformationToFrame(const Vec3 &ray_vel,
                                                const Vec3 &rf_vel);
  // Add Lorentz transformation from a frame to transform.
  // Args:
  //  ray_vel: velocity of light ray whose color spectrum is to be transformed
  //  rf_vel: velocity of reference frame from which we transform
  SpectrumTransform &ApplyTransformationFromFrame(const Vec3 &ray_vel,
                                                  const Vec3 &rf_vel);

 private:
  // overall brightness multiplier
  float brightness_;
  // rescale factor for wavelengths in spectrum
  float wavelength_rescale_;
  // absorption spectrum to be applied - given in camera frame, treated as
  // fraction of light to be absorbed as a function of wavelength
  std::vector<Spectrum> camera_frame_absorptions_;
};

#endif
