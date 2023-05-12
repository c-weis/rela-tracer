// Copyright 2023 Christoph Weis
#include <cstdio>
#include <iostream>

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

// Helpful test vars
Vec3 x(1, 0, 0);
Vec3 y(0, 1, 0);
Vec3 z(0, 0, 1);

Scene ball_in_box(int width, int height) {
  /*
    BALL IN BOX
    */
  Scene scene;

  // Camera
  Vec3 cam_pos = kZero3;
  Vec3 forward = -z;
  Vec3 right = x;
  Vec3 up = Cross(right, forward);
  Vec3 backward = -forward;
  Vec3 left = -right;
  Vec3 down = -up;

  Camera cam(cam_pos, forward, right, width, height);
  scene.AddCamera(cam);

  // Materials
  // Material bright_red_light = kBrightRed;
  // Component bright_green_light = kBrightGreen;
  // Component bright_blue_light = kBrightBlue;
  // Component dim_blue_light = kDimBlue;
  // Component lambertian = Lambertian(0.7f);
  // Component mirror = kMirror;
  // Component glass = kGlass;
  Material glass = kGlass;
  Material whiteLight = Metal().SetProperty("emission", kWhite);

  // Auxiliary screen & direction variables
  float aspect = static_cast<float>(width) / height;
  Vec3 screen_top = (right ^ forward) / aspect;

  // Add objects
  Vec3 sphere_vel(0, 0, 0);
  Vec3 sphere_O = 2 * forward;
  float sphere_rad = 1.0f;
  Sphere sphere(glass, sphere_O, sphere_rad, 0, sphere_vel);
  scene.Add(sphere);

  /* COMMENTED OUT TO ALLOW COMPILATION BEFORE IMPLEMENTATION OF COLORS
    Vec3 light_sphere_O = cam_pos + 0.5 * right + 1.5f * screen_top;
    float light_sphere_rad = 0.6f;
    Sphere light_sphere(bright_blue_light, light_sphere_O, light_sphere_rad);
    // scene.Add(light_sphere);
    //  Sphere small_mirror_sphere(&mirror, light_sphere_O, light_sphere_rad);
    //  scene.Add(small_mirror_sphere);

    Vec3 FDL = sphere_O - 3 * (forward + up + right);
    Vec3 BDL = FDL + 6 * forward;
    Vec3 FUL = FDL + 6 * up;
    Vec3 FDR = FDL + 6 * right;
    Box surrounding_box(lambertian, FDL, BDL, FUL, FDR);
    surrounding_box.back.SetMaterial(mirror);
    surrounding_box.down.SetMaterial(mirror);
    surrounding_box.right.SetMaterial(glass);
    surrounding_box.up.SetMaterial(bright_blue_light);
    surrounding_box.left.SetMaterial(bright_red_light);
    surrounding_box.front.SetMaterial(bright_green_light);
    scene.Add(surrounding_box);
  */

  return scene;
}

Scene moving_sphere(int width, int height, float speed) {
  /*
    MOVING SPHERE
  */
  Scene scene;

  // Camera
  Vec3 cam_pos = kZero3;
  Vec3 forward = -z;
  Vec3 right = x;
  Vec3 up = Cross(right, forward);
  Vec3 backward = -forward;
  Vec3 left = -right;
  Vec3 down = -up;

  Camera cam(cam_pos, forward, right, width, height);
  scene.AddCamera(cam);

  Vec3 skylight_direction = down.NormalizedNonzero();
  scene.SetSkylight(kLightYellow, skylight_direction)
      .SetAmbientBackground(kDimCyan);

  // Materials
  Material fuzzy_metal = Metal(0.8f, kBlack, 0.4f);
  Material mirror = Mirror();
  Material white_light = Metal().SetProperty("emission", kWhite);
  Material green_glowing_mirror =
      Metal().SetProperty("emission", Spectrum({GaussianData(0.5f, 530, 20)}));

  // Add objects
  // Vec3 sphere_vel = 0.9 * right;
  Vec3 sphere_vel = speed * right;
  Vec3 sphere_O = 3 * forward;
  float sphere_rad = 1.0f;
  Sphere sphere(green_glowing_mirror, sphere_O, sphere_rad, 0, sphere_vel);
  scene.Add(sphere);

  Vec3 plane_O = sphere_O + down * 3 * sphere_rad;
  Vec3 plane_A = plane_O + right;
  Vec3 plane_B = plane_O + forward;
  Plane plane(fuzzy_metal, plane_O, plane_A, plane_B);
  scene.Add(plane);

  float mirror_side = 8.0f;
  Vec3 mirror_O = 10 * forward - mirror_side / 2.0f * (up + right + backward);
  Vec3 mirror_A = mirror_O + mirror_side * (right + backward);
  Vec3 mirror_B = mirror_O + mirror_side * up;
  Parallelogram mirror_plane(mirror, mirror_O, mirror_A, mirror_B);
  scene.Add(mirror_plane);

  return scene;
}

Scene glass_ball(int width, int height) {
  Scene scene;

  Vec3 forward = y;
  Vec3 right = x;
  Vec3 up = right ^ forward;
  Camera cam(kZero3 - forward + up, forward, right, width, height);
  scene.AddCamera(cam);

  Material lambertian = Material(0.9f, 0.7f);
  Material clean_metal = Metal(0.7f);
  Material glass = kGlass;

  //scene.SetSkylight(kBrightYellow, -up);

  Vec3 O1 = kZero3;
  scene.Add(Plane(clean_metal, O1, O1 + right, O1 + forward));

  float r2 = 0.7f;
  Vec3 O2 = O1 + 2 * forward + r2 * up;
  scene.Add(Sphere(glass, O2, r2));

  float r3 = 0.5f;
  Vec3 O3 = O2 + (r2 + r3) * right - r2 * up + r3 * up + r3 / 2.0f * forward;
  scene.Add(Sphere(Metal().SetProperty("emission", kBrightGreen), O3, r3));

  float r4 = 0.5f;
  Vec3 O4 = O2 - (r2 + r4) * right - r2 * up + r4 * up;
  scene.Add(Sphere(Metal().SetProperty("emission", kBrightBlue), O4, r4));

  return scene;
}

int main(int argc, char *argv[]) {
  // Basic settings
  const int width = 500;
  const int height = 500;
  const int rays_per_pixel = 100;
  const int depth = 8;

  // Scene scene = moving_sphere(width, height, 0.95);
  Scene scene = glass_ball(width, height);

  // Raytracer
  Tracer tracey(scene);
  tracey.RenderImage(rays_per_pixel, depth, 0, 5, "glass_ball");
  // tracey.RenderFilm(rays_per_pixel, depth, 0, 25, "glass_ball", 0.0f, 1.0f,
  // 0.01f);

  return EXIT_SUCCESS;
}

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
