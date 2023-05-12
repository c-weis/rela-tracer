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

// Color constants
const Spectrum kBrightRed = BrightNarrowBand(680);
const Spectrum kBrightOrange = BrightNarrowBand(605);
const Spectrum kBrightYellow = BrightNarrowBand(570);
const Spectrum kBrightGreen = BrightNarrowBand(530);
const Spectrum kBrightCyan = BrightNarrowBand(490);
const Spectrum kBrightBlue = BrightNarrowBand(460);
const Spectrum kBrightPurple = BrightNarrowBand(410);
const Spectrum kLightRed = LightNarrowBand(680);
const Spectrum kLightOrange = LightNarrowBand(605);
const Spectrum kLightYellow = LightNarrowBand(580);
const Spectrum kLightGreen = LightNarrowBand(530);
const Spectrum kLightCyan = LightNarrowBand(490);
const Spectrum kLightBlue = LightNarrowBand(460);
const Spectrum kLightPurple = LightNarrowBand(410);
const Spectrum kDimRed = DimNarrowBand(680);
const Spectrum kDimOrange = DimNarrowBand(605);
const Spectrum kDimYellow = DimNarrowBand(580);
const Spectrum kDimGreen = DimNarrowBand(530);
const Spectrum kDimCyan = DimNarrowBand(490);
const Spectrum kDimBlue = DimNarrowBand(460);
const Spectrum kDimPurple = DimNarrowBand(410);

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
  float std_dev_4 = sqrtf(second_moment_4 - first_moment_4_normSq);

  std::cout << "Generated " << kSampleSize << " random unit vectors."
            << std::endl
            << "Mean: (" << first_moment_4.x << ", " << first_moment_4.y << ")"
            << std::endl
            << "Standard deviation: " << std_dev_4 << std::endl;

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

  Spectrum super_dim_white_emission(kWhite * 0.01);

  std::string filename1 = "images/lambda_brightness.bmp";
  std::string filename2 = "images/lambda_bandwidth.bmp";
  std::string filename3 = "images/lambda_lambda.bmp";
  std::string filename4 = "images/emit_absorb.bmp";
  std::string filename5 = "images/absorb_absorb_parallel.bmp";
  std::string filename6 = "images/absorb_absorb_series.bmp";

  SpectrumTransform no_transform;

  std::cout << "Testing narrowband emission. Varying lambda and brightness."
            << std::endl;

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
            GaussianData(narrowband_amplitude, emit_lambda, narrowband_sig) *
            brightness);
        RGBData pixel_RGB =
            no_transform.ColorFrom(emission).ToRGB(rescale_factor);
        test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
      }
    }
  }

  std::cout << "Outputting to " << filename1 << "." << std::endl;
  test_image.save_image(filename1);
  std::cout << "----" << std::endl;

  std::cout
      << "Testing generic gaussian emission. Varying mean lambda and bandwidth."
      << std::endl;

  test_image = BitmapImage(lambda_side, bandwidth_side);

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
        Spectrum emission(GaussianData(1.0f / bandwidth, lambda, bandwidth));
        RGBData pixel_RGB =
            no_transform.ColorFrom(emission).ToRGB(rescale_factor);
        test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
      }
    }
  }

  std::cout << "Outputting to " << filename2 << "." << std::endl;
  test_image.save_image(filename2);
  std::cout << "----" << std::endl;

  std::cout
      << "Testing color mixing of narrow band emissions. Varying lambda_x "
         "and lambda_y."
      << std::endl;

  test_image = BitmapImage(lambda_side, lambda_side);

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
            {GaussianData(narrowband_amplitude / 2, lambda_x, narrowband_sig),
             GaussianData(narrowband_amplitude / 2, lambda_y, narrowband_sig)});
        RGBData pixel_RGB =
            no_transform.ColorFrom(emission).ToRGB(rescale_factor);
        test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
      }
    }
  }

  std::cout << "Outputting to " << filename3 << "." << std::endl;
  test_image.save_image(filename3);
  std::cout << "----" << std::endl;

  std::cout << "Testing narrow band emission and absorption (with higher "
               "bandwidth). Varying lambda_emit "
               "and lambda_absorb."
            << std::endl;

  test_image = BitmapImage(lambda_side, lambda_side);

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
            GaussianData(narrowband_amplitude, lambda_emit, narrowband_sig));
        Spectrum absorption(GaussianData(narrowband_amplitude, lambda_absorb,
                                         absorptionband_sig));
        SpectrumTransform new_transform(no_transform);
        new_transform.ApplyAbsorption(absorption);
        RGBData pixel_RGB =
            new_transform.ColorFrom(emission).ToRGB(rescale_factor);
        test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
      }
    }
  }

  std::cout << "Outputting to " << filename4 << "." << std::endl;
  test_image.save_image(filename4);
  std::cout << "----" << std::endl;

  std::cout << "Testing single absorption with multiple modes. Varying "
               "lambda_absorb_x "
               "and lambda_absorb_y."
            << std::endl;

  test_image = BitmapImage(lambda_side, lambda_side);

  for (int y = 0; y < lambda_side; y++) {
    float lambda_absorb_y = min_lambda + d_lambda * y;
    for (int x = 0; x < lambda_side; x++) {
      float lambda_absorb_x = min_lambda + d_lambda * x;

      if ((static_cast<int>(floorf(lambda_absorb_x)) % 100 < d_lambda) ||
          (static_cast<int>(floorf(lambda_absorb_y)) % 100 < d_lambda)) {
        // place white tick
        test_image.set_pixel(x, y, 0, 0, 0);
      } else {
        // place actual color: emission at lambda_emit, absorption at
        // lambda_absorb
        Spectrum absorption({GaussianData(narrowband_amplitude, lambda_absorb_x,
                                          absorptionband_sig),
                             GaussianData(narrowband_amplitude, lambda_absorb_y,
                                          absorptionband_sig)});
        SpectrumTransform new_transform(no_transform);
        new_transform.ApplyAbsorption(absorption);
        RGBData pixel_RGB = new_transform.ColorFrom(super_dim_white_emission)
                                .ToRGB(rescale_factor);
        test_image.set_pixel(x, y, pixel_RGB.R, pixel_RGB.G, pixel_RGB.B);
      }
    }
  }

  std::cout << "Outputting to " << filename5 << "." << std::endl;
  test_image.save_image(filename5);
  std::cout << "----" << std::endl;

  std::cout << "Testing two absorptions with a single mode each. Varying "
               "lambda_absorb_x "
               "and lambda_absorb_y."
            << std::endl;

  test_image = BitmapImage(lambda_side, lambda_side);

  for (int y = 0; y < lambda_side; y++) {
    float lambda_absorb_y = min_lambda + d_lambda * y;
    for (int x = 0; x < lambda_side; x++) {
      float lambda_absorb_x = min_lambda + d_lambda * x;

      if ((static_cast<int>(floorf(lambda_absorb_x)) % 100 < d_lambda) ||
          (static_cast<int>(floorf(lambda_absorb_y)) % 100 < d_lambda)) {
        // place white tick
        test_image.set_pixel(x, y, 0, 0, 0);
      } else {
        // place actual color: emission at lambda_emit, absorption at
        // lambda_absorb
        Spectrum absorption1(GaussianData(narrowband_amplitude, lambda_absorb_x,
                                          absorptionband_sig));
        Spectrum absorption2(GaussianData(narrowband_amplitude, lambda_absorb_y,
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

  std::cout << "Outputting to " << filename6 << "." << std::endl;
  test_image.save_image(filename6);
  std::cout << "----" << std::endl;

  return true;
}

bool Tester::TestLorentzTransforms() {}

bool Tester::TestIntersections() {}

bool Tester::RunAllTests() { return TestRandomness() && TestColors(); }

/* PASSED: Test Lambertian scattering
// (avg should give -2/3 along chosen norm, avg_normsq should be 1)
std::vector<Vec3> l_scatters;
const int nr_scatters = 10000;
l_scatters.reserve(nr_scatters);
std::mt19937 r_gen(42);
for (int i = 0; i < nr_scatters; i++) {
  l_scatters.push_back(
      LambertianScatter::LambertianInverseRay(Vec3(1, 0, 0), &r_gen));
}
Vec3 avg = std::reduce(l_scatters.cbegin(), l_scatters.cend(), kZero3) /
           static_cast<float>(nr_scatters);
float avg_normsq = 0;
std::for_each(l_scatters.cbegin(), l_scatters.cend(),
              [&avg_normsq](Vec3 vec) { avg_normsq += vec.NormSq(); });
avg_normsq /= nr_scatters;
std::cout << "Lambert scatter: "
          << "avg: " << avg << " "
          << "normsq: " << avg_normsq << " "
          << std::endl;
*/

/* PASSED: Line transformation test
Line testline1 = Line(Vec4(0,1,2,3), Vec3(-0.5,0.5,0.5));
Line testline2 = Line(Vec4(-1,11,-2,5), Vec3(0.5,-0.5,0.5));
std::cout << "Line 1 in frame 0: " << testline1 << std::endl;
std::cout << "Line 1 in frame 1: " << testline1.in_frame(testline1) <<
std::endl; std::cout << "Line 1 in frame 0: " <<
testline1.in_frame(testline1).TransformedFromFrame(testline1) << std::endl;
std::cout << "Line 1 in frame 2: " << testline1.in_frame(testline2) <<
std::endl; std::cout
<< "Line 2 in frame 0: " << testline2 << std::endl; std::cout << "Line 2 in
frame 1: " << testline2.in_frame(testline1) << std::endl; std::cout << "Line 2
in frame 0: " << testline2.in_frame(testline1).TransformedFromFrame(testline1)
<< std::endl; std::cout << "Line 2 in frame 2: " <<
testline2.in_frame(testline2)
<< std::endl;
*/

/* PASSED: Image Line inspection
LineList img_rays = cam.ImageRays();
std::cout << "Number of rays: " << img_rays.size() << std::endl;
for (Line r : img_rays) {
  std::printf(
      "origin: (%1.1f, %1.1f, %1.1f, %1.1f) vel: (%1.1f, %1.1f, %1.1f)\n",
      r.origin.t, r.origin.r.x, r.origin.r.y, r.origin.r.z,
      r.vel.x, r.vel.y, r.vel.z);
}
*/

/* PASSED: Basic sphere test
std::cout << "Sphere intersection test: " << std::endl;
LineList test_rays{Line(kZero4, xvel),
                  Line(Vec4(0, 0, 0.9, 0), xvel),
                  Line(Vec4(0, 0, 0, 0.9), xvel),
                  Line(Vec4(0, 0, -0.9, 0), xvel),
                  Line(Vec4(0, 0, 0, -0.9), xvel),
                  Line(Vec4(0, 0, 1.1, 0), xvel),
                  Line(Vec4(0, 0, 0, 1.1), xvel),
                  Line(Vec4(0, 0, -1.1, 0), xvel),
                  Line(Vec4(0, 0, 0, -1.1), xvel),
                  };

for (Line ray : test_rays) {
  printf("Test Line (y,z) = (%1.1f, %1.1f) in x-direction ",
         ray.origin.r.y, ray.origin.r.z);
  if (s.intersect(ray)) {
    printf("hits.");
  } else {
    printf("doesn't hit.");
  }
  printf("\n");
}

*/

int main() {
  Tester test;

  return test.RunAllTests();
}
