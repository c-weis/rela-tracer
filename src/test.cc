// Copyright 2023 Christoph Weis
#include "include/test.h"

#include <cstdio>
#include <iostream>

#include "include/bitmap_image.h"
#include "include/colors.h"
#include "include/math.h"
#include "include/tracer.h"

// TODO(c): figure out how to clean this up:
// include these in the lib code somehow

bool Tester::TestRandomness() {
  const int kSampleSize = 10000;

  std::cout << "Testing randomness functions." << std::endl;

  std::cout << "---" << std::endl;
  std::cout << "Generating uniform random reals in [0, 1)." << std::endl;

  float first_moment_1 = 0;
  float second_moment_1 = 0;
  for (int i = 0; i < kSampleSize; i++) {
    float real = RandomReal();
    if (real < 0 || real >= 1) {
      std::cout << "Error: produced random value " << real << " outside range."
                << std::endl;
      return false;
    }
    first_moment_1 += real;
    second_moment_1 += real * real;
  }
  first_moment_1 /= kSampleSize;
  second_moment_1 /= kSampleSize;
  float std_dev_1 = sqrtf(second_moment_1 - first_moment_1 * first_moment_1);

  std::cout << "Generated " << kSampleSize << " random values." << std::endl
            << "Mean: " << first_moment_1 << std::endl
            << "Standard deviation: " << std_dev_1 << std::endl;
  std::cout << "---" << std::endl;

  std::cout << "Generating random vectors in unit ball." << std::endl;
  Vec3 first_moment_2(0, 0, 0);
  float second_moment_2 = 0;
  for (int i = 0; i < kSampleSize; i++) {
    Vec3 random_vec = RandomVectorInUnitBall();
    if (random_vec.NormSq() >= 1 || random_vec.NormSq() == 0) {
      std::cout << "Error: produced random vector " << random_vec
                << " with norm " << random_vec.Norm() << " outside range."
                << std::endl;
      return false;
    }
    first_moment_2 = first_moment_2 + random_vec;
    second_moment_2 += random_vec.NormSq();
  }
  first_moment_2 = first_moment_2 / kSampleSize;
  second_moment_2 = second_moment_2 / kSampleSize;
  float std_dev_2 =
      sqrtf(second_moment_2 - Dot3(first_moment_2, first_moment_2));

  std::cout << "Generated " << kSampleSize << " random vectors." << std::endl
            << "Mean: " << first_moment_2 << std::endl
            << "Standard deviation: " << std_dev_2 << std::endl;
  std::cout << "---" << std::endl;

  std::cout << "Generating random unit vectors." << std::endl;
  float kEpsilon = 1e-6;
  Vec3 first_moment_3(0, 0, 0);
  float second_moment_3 = 0;
  for (int i = 0; i < kSampleSize; i++) {
    Vec3 random_vec = RandomUnitVector();
    if (random_vec.NormSq() > 1 + kEpsilon ||
        random_vec.NormSq() < 1 - kEpsilon) {
      std::cout << "Error: produced random vector " << random_vec
                << " with norm " << random_vec.Norm()
                << " outside allowable range." << std::endl;
      return false;
    }
    first_moment_3 = first_moment_3 + random_vec;
    second_moment_3 += random_vec.NormSq();
  }
  first_moment_3 = first_moment_3 / kSampleSize;
  second_moment_3 = second_moment_3 / kSampleSize;
  float std_dev_3 =
      sqrtf(second_moment_3 - Dot3(first_moment_3, first_moment_3));

  std::cout << "Generated " << kSampleSize << " random unit vectors."
            << std::endl
            << "Mean: " << first_moment_3 << std::endl
            << "Standard deviation: " << std_dev_3 << std::endl;
  std::cout << "---" << std::endl;

  std::cout << "Generating random vectors in unit disk." << std::endl;
  Vec2 first_moment_4(0, 0);
  float second_moment_4 = 0;
  for (int i = 0; i < kSampleSize; i++) {
    Vec2 random_vec = RandomVectorInUnitDisk();
    float normSq = random_vec.x * random_vec.x + random_vec.y * random_vec.y;
    if (normSq >= 1) {
      std::cout << "Error: produced random vector (" << random_vec.x << ", "
                << random_vec.y << ") with norm " << sqrtf(normSq)
                << " outside range." << std::endl;
      return false;
    }
    first_moment_4.x = first_moment_4.x + random_vec.x;
    first_moment_4.y = first_moment_4.y + random_vec.y;
    second_moment_4 += normSq;
  }
  first_moment_4.x = first_moment_4.x / kSampleSize;
  first_moment_4.y = first_moment_4.y / kSampleSize;
  float first_moment_4_normSq =
      first_moment_4.x * first_moment_4.x + first_moment_4.y * first_moment_4.y;
  second_moment_4 = second_moment_4 / kSampleSize;

  std::cout << "Generated " << kSampleSize << " random unit vectors."
            << std::endl
            << "Mean: (" << first_moment_4.x << ", " << first_moment_4.y << ")"
            << std::endl
            << "Mean norm squared: " << second_moment_4 << std::endl;

  std::cout << "---" << std::endl;
  std::cout << "Tests completed successfully." << std::endl;

  return true;
}

bool Tester::TestColors() {
  const float min_lambda = 300;  // nm
  const float max_lambda = 800;  // nm
  const float d_lambda = 0.5;
  const unsigned int lambda_side = (max_lambda - min_lambda) / d_lambda;

  const float narrowband_sig = 0.5f;        // nm
  const float narrowband_amplitude = 1.0f;  // nm
  const float absorptionband_sig = 50.0f;

  const float min_brightness = 0.0f;
  const float max_brightness = 1.0f;
  const float d_brightness = 0.001f;
  const float rescale_factor = 255.0f;
  const unsigned int brightness_side =
      (max_brightness - min_brightness) / d_brightness;

  const float min_bandwidth = 0.0f;
  const float max_bandwidth = 300.0f;
  const float d_bandwidth = 0.5f;
  const unsigned int bandwidth_side =
      (max_bandwidth - min_bandwidth) / d_bandwidth;

  Spectrum super_dim_white_emission(Spectrum::White * 0.01);

  std::string folder = "images/test/";

  SpectrumTransform no_transform;

  // Part 1: emission - vary lambda and brightness
  {
    std::cout << "Testing narrowband emission. Varying lambda and brightness."
              << std::endl;

    std::string filename(folder + "lambda_brightness.bmp");

    BitmapImage test_image(lambda_side, brightness_side);

    for (int y = 0; y < brightness_side; y++) {
      float brightness = min_brightness + d_brightness * y;
      for (int x = 0; x < lambda_side; x++) {
        float emit_lambda = min_lambda + d_lambda * x;

        if (static_cast<int>(floorf(emit_lambda)) % 100 < d_lambda) {
          // place white tick
          test_image.set_pixel(x, y, 255, 255, 255);
        } else if (static_cast<int>(floorf(emit_lambda)) % 50 < d_lambda) {
          // place grey tick
          test_image.set_pixel(x, y, 100, 100, 100);
        } else {
          // place actual color
          Spectrum emission(
              Gaussian(narrowband_amplitude, emit_lambda, narrowband_sig) *
              brightness);
          RGBData pixel_RGB =
              no_transform.ColorFrom(emission).ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
    std::cout << "----" << std::endl;
  }

  // Part 2: emission - vary lambda and bandwidth
  {
    std::cout << "Testing generic gaussian emission. Varying mean lambda and "
                 "bandwidth."
              << std::endl;

    std::string filename(folder + "lambda_bandwidth.bmp");

    BitmapImage test_image(lambda_side, bandwidth_side);

    for (int y = 0; y < bandwidth_side; y++) {
      float bandwidth = min_bandwidth + d_bandwidth * y;
      for (int x = 0; x < lambda_side; x++) {
        float lambda = min_lambda + d_lambda * x;

        if (static_cast<int>(floorf(lambda)) % 100 < d_lambda) {
          // place white tick
          test_image.set_pixel(x, y, 255, 255, 255);
        } else if (static_cast<int>(floorf(lambda)) % 50 < d_lambda) {
          // place grey tick
          test_image.set_pixel(x, y, 100, 100, 100);
        } else {
          // place actual color
          Spectrum emission(Gaussian(1.0f / bandwidth, lambda, bandwidth));
          RGBData pixel_RGB =
              no_transform.ColorFrom(emission).ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
    std::cout << "----" << std::endl;
  }

  // Part 3: emission color mixing - vary lambda and lambda
  {
    std::cout
        << "Testing color mixing of narrow band emissions. Varying lambda_x "
           "and lambda_y."
        << std::endl;

    std::string filename(folder + "lambda_lambda.bmp");

    BitmapImage test_image(lambda_side, lambda_side);

    for (int y = 0; y < lambda_side; y++) {
      float lambda_y = min_lambda + d_lambda * y;
      for (int x = 0; x < lambda_side; x++) {
        float lambda_x = min_lambda + d_lambda * x;

        if ((static_cast<int>(floorf(lambda_x)) % 100 < d_lambda) ||
            (static_cast<int>(floorf(lambda_y)) % 100 < d_lambda)) {
          // place white tick
          test_image.set_pixel(x, y, 255, 255, 255);
        } else {
          // place actual color: a mix of lambda_x and lambda_y narrowbands
          Spectrum emission(
              {Gaussian(narrowband_amplitude / 2, lambda_x, narrowband_sig),
               Gaussian(narrowband_amplitude / 2, lambda_y, narrowband_sig)});
          RGBData pixel_RGB =
              no_transform.ColorFrom(emission).ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
    std::cout << "----" << std::endl;
  }

  // Part 4: emission-absorption: vary lambda_emit & lambda_absorb
  {
    std::cout << "Testing narrow band emission and absorption (with higher "
                 "bandwidth). Varying lambda_emit "
                 "and lambda_absorb."
              << std::endl;

    std::string filename(folder + "emit_absorb.bmp");

    BitmapImage test_image(lambda_side, lambda_side);

    for (int y = 0; y < lambda_side; y++) {
      float lambda_absorb = min_lambda + d_lambda * y;
      for (int x = 0; x < lambda_side; x++) {
        float lambda_emit = min_lambda + d_lambda * x;

        if ((static_cast<int>(floorf(lambda_emit)) % 100 < d_lambda) ||
            (static_cast<int>(floorf(lambda_absorb)) % 100 < d_lambda)) {
          // place white tick
          test_image.set_pixel(x, y, 255, 255, 255);
        } else {
          // place actual color: emission at lambda_emit, absorption at
          // lambda_absorb
          Spectrum emission(
              Gaussian(narrowband_amplitude, lambda_emit, narrowband_sig));
          Spectrum absorption(Gaussian(narrowband_amplitude, lambda_absorb,
                                       absorptionband_sig));
          SpectrumTransform new_transform(no_transform);
          new_transform.ApplyAbsorption(absorption);
          RGBData pixel_RGB =
              new_transform.ColorFrom(emission).ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
    std::cout << "----" << std::endl;
  }

  // Part 5: double absorption in parallel - vary lambda and lambda
  {
    std::cout << "Testing single absorption with multiple modes. Varying "
                 "lambda_absorb_x "
                 "and lambda_absorb_y."
              << std::endl;

    std::string filename(folder + "absorb_absorb_parallel.bmp");

    BitmapImage test_image(lambda_side, lambda_side);

    for (int y = 0; y < lambda_side; y++) {
      float lambda_absorb_y = min_lambda + d_lambda * y;
      for (int x = 0; x < lambda_side; x++) {
        float lambda_absorb_x = min_lambda + d_lambda * x;

        if ((static_cast<int>(floorf(lambda_absorb_x)) % 100 < d_lambda) ||
            (static_cast<int>(floorf(lambda_absorb_y)) % 100 < d_lambda)) {
          // place white tick
          test_image.set_pixel(x, y, 0, 0, 0);
        } else {
          // place actual color: dim white with absorption at
          // lambda_x and lambda_y in parallel
          Spectrum absorption({Gaussian(narrowband_amplitude, lambda_absorb_x,
                                        absorptionband_sig),
                               Gaussian(narrowband_amplitude, lambda_absorb_y,
                                        absorptionband_sig)});
          SpectrumTransform new_transform(no_transform);
          new_transform.ApplyAbsorption(absorption);
          RGBData pixel_RGB = new_transform.ColorFrom(super_dim_white_emission)
                                  .ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
    std::cout << "----" << std::endl;
  }

  // Part 6: double absorption in series: vary lambda and lambda
  {
    std::cout << "Testing two absorptions with a single mode each. Varying "
                 "lambda_absorb_x "
                 "and lambda_absorb_y."
              << std::endl;

    std::string filename(folder + "absorb_absorb_series.bmp");

    BitmapImage test_image(lambda_side, lambda_side);

    for (int y = 0; y < lambda_side; y++) {
      float lambda_absorb_y = min_lambda + d_lambda * y;
      for (int x = 0; x < lambda_side; x++) {
        float lambda_absorb_x = min_lambda + d_lambda * x;

        if ((static_cast<int>(floorf(lambda_absorb_x)) % 100 < d_lambda) ||
            (static_cast<int>(floorf(lambda_absorb_y)) % 100 < d_lambda)) {
          // place white tick
          test_image.set_pixel(x, y, 0, 0, 0);
        } else {
          // place actual color: dim white with absorption at
          // lambda_x and lambda_y in series
          Spectrum absorption1(Gaussian(narrowband_amplitude, lambda_absorb_x,
                                        absorptionband_sig));
          Spectrum absorption2(Gaussian(narrowband_amplitude, lambda_absorb_y,
                                        absorptionband_sig));
          SpectrumTransform new_transform(no_transform);
          new_transform.ApplyAbsorption(absorption1);
          new_transform.ApplyAbsorption(absorption2);
          RGBData pixel_RGB = new_transform.ColorFrom(super_dim_white_emission)
                                  .ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
    std::cout << "----" << std::endl;
  }

  return true;
}

float kMaxEntry = 1e2;
Vec4 RandomPos4() {
  return kMaxEntry * Vec4(2 * RandomReal() - 1, 2 * RandomReal() - 1,
                          2 * RandomReal() - 1, 2 * RandomReal() - 1);
}

Vec3 RandomVel3() {
  // Random vector in unit ball, heavily biased towards center
  return RandomVectorInUnitBall() * RandomReal() * RandomReal();
}

float EuclNormSq4(Vec4 vec) { return vec.SpaceNormSq() + vec.t * vec.t; }

bool Tester::TestLineLorentzTransforms() {
  float kTestEpsilon = 1e-6;

  std::cout << "Testing 4-line transformations: " << std::endl
            << "----" << std::endl;

  // Part 1: to-and-from
  {
    std::cout << "Transform to-and-from:" << std::endl;
    Vec4 o = RandomPos4();
    Vec3 v = RandomVel3();

    Vec4 rf_o = RandomPos4();
    Vec3 rf_v = RandomVel3();

    Line l(o, v);
    Line rf(rf_o, rf_v);
    std::cout << "Random line:" << l << std::endl;
    std::cout << "Random reference frame worldline rf: " << rf << std::endl;

    Line l_l = l.TransformedToFrame(l);
    Line l_rf = l.TransformedToFrame(rf);
    Line l_rf_0 = l_rf.TransformedFromFrame(rf);
    std::cout << "Line in its own rest frame: " << l_l << std::endl
              << "Line in reference frame rf: " << l_rf << std::endl
              << "Line transformed back to standard frame: " << l_rf_0
              << std::endl;

    // A line should be at rest at origin in its rest frame.
    float eucl_normSq = EuclNormSq4(l_l.origin);
    float vel_normSq = l_l.vel.NormSq();
    if (eucl_normSq > kTestEpsilon) {
      std::cout << "Line not at origin in its rest frame!" << std::endl;
      return false;
    }
    if (vel_normSq > kTestEpsilon) {
      std::cout << "Line not at rest in its rest frame!" << std::endl;
      return false;
    }
    if (EuclNormSq4(l.origin - l_rf_0.origin) >
            kTestEpsilon * EuclNormSq4(l.origin) ||
        (l.vel - l_rf_0.vel).NormSq() > kTestEpsilon * l.vel.NormSq()) {
      std::cout << "Line does not return to original value under to-and-from "
                   "transformation."
                << std::endl;
      return false;
    }
    std::cout << "----" << std::endl;
  }

  // Part 2: linearity
  {
    std::cout << "Check linearity:" << std::endl;

    Vec3 rf_vel = RandomVel3();

    Vec4 pos1 = RandomPos4();
    Vec4 delta = RandomPos4();
    Vec4 pos2 = pos1 + delta;

    std::cout << "Random positions: " << std::endl
              << "pos1: " << pos1 << std::endl
              << "pos2 = pos1 + delta: " << pos2 << std::endl;
    std::cout << "Random rf velocity: " << rf_vel << std::endl;

    Vec4 pos1_rf = pos1.TransformedToFrame(rf_vel);
    Vec4 delta_rf = delta.TransformedToFrame(rf_vel);
    Vec4 pos2_rf = pos2.TransformedToFrame(rf_vel);

    Vec4 sum_rf = pos1_rf + delta_rf;

    std::cout << "pos1 + delta transformed to rf: " << sum_rf << std::endl
              << "pos2 transformed to rf: " << pos2_rf << std::endl;

    if (EuclNormSq4(sum_rf - pos2_rf) > kTestEpsilon * EuclNormSq4(pos2_rf)) {
      std::cout
          << "Linearity test failed. Lorentz transformations should be linear."
          << std::endl;
    }
    std::cout << "----" << std::endl;
  }

  // Part 3: transitivity
  {
    std::cout << "Check transform transitivity:" << std::endl;
    // keep speed low to avoid use of large numbers
    float max_speed = 0.3;
    Vec4 o = RandomPos4();
    Vec3 v = RandomVel3() * max_speed;
    Vec4 o1 = RandomPos4();
    Vec3 v1 = RandomVel3() * max_speed;
    Vec4 o2 = RandomPos4();
    Vec3 v2 = RandomVel3() * max_speed;

    Line l(o, v);
    Line rf1(o1, v1);
    Line rf2(o2, v2);

    std::cout << "Random line: " << l << std::endl;
    std::cout << "Random rf worldlines: " << std::endl
              << "rf1: " << rf1 << std::endl
              << "rf2: " << rf2 << std::endl;

    Line l_1 = l.TransformedToFrame(rf1);
    Line rf2_1 = rf2.TransformedToFrame(rf1);
    Line l_1_2 = l_1.TransformedToFrame(rf2_1);
    Line l_1_2_0 = l_1_2.TransformedFromFrame(rf2);

    std::cout << "Line transformed to rf1: " << l_1 << std::endl
              << "rf2 transformed to rf1: " << rf2_1 << std::endl
              << "Line transformed to rf2 via rf1: " << l_1_2 << std::endl
              << "Line transformed back to standard frame: " << l_1_2_0
              << std::endl;

    if (EuclNormSq4(l_1_2_0.origin - l.origin) >
            kTestEpsilon * EuclNormSq4(l.origin) ||
        (l_1_2_0.vel - l.vel).NormSq() > kTestEpsilon * l.vel.NormSq()) {
      std::cout
          << "Transitivity failed: line not preserved under round-trip."
          << std::endl
          << "This is likely a result of limited precision. Continuing tests."
          << std::endl;
    }
    std::cout << "----" << std::endl;
  }

  return true;
}

bool Tester::TestColorLorentzTransforms() {
  const float min_lambda = 0;     // nm
  const float max_lambda = 1000;  // nm
  const float d_lambda = 1;
  const unsigned int lambda_side = (max_lambda - min_lambda) / d_lambda;

  const float standard_lambda = (max_lambda + min_lambda) / 2.0f;  // nm

  const float narrowband_sig = 0.5f;        // nm
  const float narrowband_amplitude = 1.0f;  // nm
  const float absorptionband_sig = 50.0f;

  const float rescale_factor = 255.0f;

  const float min_speed = -1.0f;
  const float max_speed = 1.0f;
  const float d_speed = 4e-3f;
  const unsigned int speed_side = (max_speed - min_speed) / d_speed;

  const float min_angle = -kPi;
  const float max_angle = kPi;
  const float d_angle = kPi / 800.0f;
  const unsigned int angle_side = (max_angle - min_angle) / d_angle;

  Spectrum super_dim_white_emission(Spectrum::White * 0.01);
  Vec3 ray_vel(-1, 0, 0);
  SpectrumTransform no_transform;

  std::string folder = "images/test/";

  std::cout << "Testing light transforms: " << std::endl;
  std::cout << "----" << std::endl;

  // Part 1: A moving emitter - vary wavelength and speed
  {
    std::cout << "Emission from a moving body. Varying wavelength and relative "
                 "speed: "
              << std::endl;

    std::string filename(folder + "lambda_speed.bmp");

    BitmapImage test_image(lambda_side, speed_side);

    for (int y = 0; y < speed_side; y++) {
      float speed = min_speed + d_speed * y;
      Vec3 frame_vel = ray_vel * speed;
      for (int x = 0; x < lambda_side; x++) {
        float lambda = min_lambda + d_lambda * x;

        if (static_cast<int>(floorf(lambda)) % 100 < d_lambda) {
          // place white tick
          test_image.set_pixel(x, y, 255, 255, 255);
        } else {
          // place actual color:
          // lambda, sent from source with relative speed
          Spectrum emission(
              Gaussian(narrowband_amplitude, lambda, narrowband_sig));
          SpectrumTransform new_transform(no_transform);
          new_transform.ApplyTransformationFromFrame(ray_vel, frame_vel);
          RGBData pixel_RGB =
              new_transform.ColorFrom(emission).ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
    std::cout << "----" << std::endl;
  }

  // Part 2: moving emitter of fixed wavelength - varying relative angle and
  // speed
  {
    std::cout << "Varying relative angle and speed: " << std::endl;
    std::cout << "Emitting narrow band beam with wavelength of "
              << standard_lambda << " nm" << std::endl;

    std::string filename(folder + "angle_speed.bmp");

    BitmapImage test_image(angle_side, speed_side);

    for (int y = 0; y < speed_side; y++) {
      float speed = min_speed + d_speed * y;
      for (int x = 0; x < angle_side; x++) {
        float angle = min_angle + d_angle * x;
        Vec3 frame_vel = speed * Vec3(cosf(angle), sinf(angle), 0);

        // place actual color:
        // standard_lambda, sent from source with specified relative speed and
        // angle
        Spectrum emission(
            Gaussian(narrowband_amplitude, standard_lambda, narrowband_sig));
        SpectrumTransform new_transform(no_transform);
        new_transform.ApplyTransformationFromFrame(ray_vel, frame_vel);
        RGBData pixel_RGB =
            new_transform.ColorFrom(emission).ToRGB(rescale_factor);
        test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);

    std::cout << "----" << std::endl;
  }

  // Part 3: [emission + absorption] in moving frame - varying absorption
  // wavelength and speed
  {
    std::cout << "Varying absorption wavelength and speed, "
              << "assuming emission happens in same frame as absorption:"
              << std::endl;

    std::string filename(folder + "absorb_speed_abframe.bmp");

    BitmapImage test_image(lambda_side, speed_side);

    for (int y = 0; y < speed_side; y++) {
      float speed = min_speed + d_speed * y;
      Vec3 frame_vel = ray_vel * speed;
      for (int x = 0; x < angle_side; x++) {
        float lambda = min_lambda + d_lambda * x;

        if (static_cast<int>(floorf(lambda)) % 100 < d_lambda) {
          // place white tick
          test_image.set_pixel(x, y, 0, 0, 0);
        } else {
          // place actual color:
          // super dim white experiencing absorption at specified relative speed
          Spectrum absorption(
              Gaussian(narrowband_amplitude, lambda, absorptionband_sig));
          SpectrumTransform new_transform(no_transform);
          new_transform.ApplyTransformationFromFrame(ray_vel, frame_vel);
          new_transform.ApplyAbsorption(absorption);

          RGBData pixel_RGB = new_transform.ColorFrom(super_dim_white_emission)
                                  .ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);

    std::cout << "----" << std::endl;
  }

  // Part 4: emission in std frame, absorption in moving frame - vary absorption
  // wavelength and speed
  {
    std::cout << "Varying absorption wavelength and speed, "
              << "assuming emission happens in standard frame:" << std::endl;

    std::string filename(folder + "absorb_speed_stdframe.bmp");

    BitmapImage test_image(lambda_side, speed_side);

    for (int y = 0; y < speed_side; y++) {
      float speed = min_speed + d_speed * y;
      Vec3 frame_vel = ray_vel * speed;
      for (int x = 0; x < angle_side; x++) {
        float lambda = min_lambda + d_lambda * x;

        if (static_cast<int>(floorf(lambda)) % 100 < d_lambda) {
          // place white tick
          test_image.set_pixel(x, y, 0, 0, 0);
        } else {
          // place actual color:
          // super dim white experiencing absorption at specified relative speed
          Spectrum absorption(
              Gaussian(narrowband_amplitude, lambda, absorptionband_sig));
          SpectrumTransform new_transform(no_transform);
          new_transform.ApplyTransformationFromFrame(ray_vel, frame_vel);
          new_transform.ApplyAbsorption(absorption);
          new_transform.ApplyTransformationToFrame(ray_vel, frame_vel);

          RGBData pixel_RGB = new_transform.ColorFrom(super_dim_white_emission)
                                  .ToRGB(rescale_factor);
          test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
        }
      }
    }

    std::cout << "Outputting to " << filename << "." << std::endl;
    test_image.save_image(filename);
  }

  return true;
}

bool Tester::TestScenes() {
  Vec3 back(1, 0, 0);
  Vec3 right(0, 1, 0);
  Vec3 down(0, 0, 1);
  Vec3 front = -back;
  Vec3 left = -right;
  Vec3 up = -down;

  int width = 500;
  int height = 500;
  int rays_per_pixel = 10;
  int depth = 3;

  Scene scene;
  scene.SetAmbientBackground(Spectrum::Cyan * 0.1f);

  Material green_glowing_bbody(0.0f, 0.0f, 0.0f, 0.0f, Spectrum::Black,
                               Spectrum::Green);
  Material lambertian(0.8f);

  Vec3 ocam = kZero3;
  Vec3 omid = ocam + back;

  float side = 0.5f;
  float diag = sqrt(3) * side;

  Camera cam(ocam, back / 2.0f, right, width, height);

  scene.Add(cam);

  float min_speed = -0.4f;
  float d_speed = 0.4f;
  int speed_steps = 3;
  float max_speed = min_speed + d_speed * speed_steps;

  std::cout << "Starting trace of test boxes: " << std::endl;
  for (int x = 0; x < speed_steps; x++) {
    float speed_right = min_speed + d_speed * x;
    for (int y = 0; y < speed_steps; y++) {
      float speed_up = min_speed + d_speed * y;
      std::cout << "speed_right: " << speed_right << ", speed_up: " << speed_up
                << ". " << std::endl;
      Scene new_scene1(scene);
      Scene new_scene2(scene);
      Box box1(omid, side / 2 * right, side / 2 * back, side / 2 * up,
               speed_right * right + speed_up * up, 0, green_glowing_bbody);
      Box box2(omid, side / 2 * right, side / 2 * back, side / 2 * up,
               speed_right * right + speed_up * up, 0, lambertian);
      new_scene1.Add(box1);
      new_scene2.Add(box2);
      Tracer testtrace1(new_scene1);
      Tracer testtrace2(new_scene2);
      testtrace1.RenderImage(
          rays_per_pixel, depth, 0, 10,
          "test/green_bbody_box_" + std::to_string(x) + "_" + std::to_string(y),
          0, 255, false);
      testtrace2.RenderImage(
          rays_per_pixel, depth, 0, 10,
          "test/lambertian_box_" + std::to_string(x) + "_" + std::to_string(y),
          0, 255, false);
    }
  }

  return true;
}

bool Tester::RunAllTests() {
  bool did_all_tests_work = TestRandomness() && TestColors() &&
                            TestLineLorentzTransforms() &&
                            TestColorLorentzTransforms() && TestScenes();
  return did_all_tests_work;
}
