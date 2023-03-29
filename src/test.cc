// Copyright 2023 Christoph Weis
#include <cstdio>
#include <iostream>

#include "include/math.h"
#include "include/scene.h"
#include "include/tracer.h"

// Helpful test vars
Vec3 x(1, 0, 0);
Vec3 y(0, 1, 0);
Vec3 z(0, 0, 1);

int main() {
  Scene scene;

  // Basic settings
  const int width = 50;
  const int height = 50;
  const int rays_per_pixel = 10;
  const int depth = 3;
  const std::vector<int> scatter_ray_counts = {5, 5, 1};

  // Camera
  Vec4 cam_pos = kZero4;
  Vec3 screen_centre = cam_pos.r + y;
  Vec3 right_edge = screen_centre + z;
  Camera cam(cam_pos, screen_centre, right_edge, width, height);
  scene.AddCamera(cam);

  // Materials
  LambertianScatter lambert;
  Material bright_red_light(nullptr, &kBrightRed);
  Material dim_green_light(nullptr, &kDimGreen);
  Material bright_blue_light(nullptr, &kBrightBlue);
  Material dim_blue_light(nullptr, &kDimBlue);
  //Material bright_red_light(&lambert, &kBrightRed);
  //Material dim_green_light(&lambert, &kDimGreen);
  //Material bright_blue_light(&lambert, &kBrightBlue);
  //Material dim_blue_light(&lambert, &kDimBlue);
  Material lambertian(&lambert, nullptr);
  Material mirror(new Mirror(), nullptr); 

  // Auxiliary screen & direction variables
  float aspect = static_cast<float>(width) / height;
  Vec3 right = right_edge - screen_centre;
  Vec3 front = screen_centre - cam_pos.r;
  Vec3 top = (right ^ front) / aspect;

  // Add objects
  Vec3 sphere_vel(0, 0, 0);
  Vec3 sphere_O = 2 * y;
  float sphere_rad = 1.0f;
  //Sphere mirror_sphere(&mirror, sphere_O, sphere_rad, 0, sphere_vel);
  //scene.Add(mirror_sphere);
  Sphere lambert_sphere(&lambertian, sphere_O, sphere_rad, 0, sphere_vel);
  scene.Add(lambert_sphere);

  Vec3 light_sphere_O = cam_pos.r + 0.5 * right + 1.5f * top;
  float light_sphere_rad = 0.3f;
  Sphere light_sphere(&bright_blue_light, light_sphere_O, light_sphere_rad);
  scene.Add(light_sphere);
  // Sphere small_mirror_sphere(&mirror, light_sphere_O, light_sphere_rad);
  // scene.Add(small_mirror_sphere);

  Vec3 box_corner_O = sphere_O - 3 * (x + y + z);
  Vec3 box_corner_A = box_corner_O + 6 * x;
  Vec3 box_corner_B = box_corner_O + 6 * y;
  Vec3 box_corner_C = box_corner_O + 6 * z;
  scene.AddBox(&lambertian, box_corner_O, box_corner_A, box_corner_B,
               box_corner_C);

  // Vec3 plane_O = sphere_O - 2 * x;
  // Plane p(&gray_light, plane_O, plane_O + x, plane_O + y);
  // scene.Add(p);

  // Line light_triangle_line(Vec4(0, screen_centre + 0.5 * right + 0.5 * top),
  // kZero3); Triangle t(&white_light, light_triangle_line, 0.5 * right, 0.5 *
  // top); scene.Add(t);

  // Raytracer
  Tracer tracey(scene);
  tracey.RenderImage(rays_per_pixel, depth, scatter_ray_counts, 0,
                     RenderMode::kOutputFile, "test");
  // tracey.RenderImage(rays_per_pixel, depth, scatter_ray_counts, 0,
  // RenderMode::kConsole);

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
