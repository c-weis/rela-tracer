// Copyright 2023 Christoph Weis
#include "include/main.h"

#include <iostream>
#include <map>
#include <vector>

#include "include/test.h"
#include "include/tracer.h"

// A helpful set of directions
const Vec3 up(0, 1, 0);
const Vec3 right(1, 0, 0);
const Vec3 forward = Cross(up, right);
const Vec3 down = -up;
const Vec3 left = -right;
const Vec3 backward = -forward;

// A static scene involving three reflective balls of differing colours.
// Scene arguments: none
Scene metal_balls(int width, int height, arglist args = {}) {
  Scene scene;

  // Materials
  Material lambertian = Material(0.9f, 0.7f);
  Spectrum absorb_red_green({Gaussian(0.5f, 550, 50), Gaussian(0.5f, 600, 50)});
  Spectrum absorb_red_blue({Gaussian(0.5f, 450, 50), Gaussian(0.5f, 600, 25)});
  Material clean_blue_metal = Material::Metal(0.7f, absorb_red_green);
  Material clean_green_metal = Material::Metal(0.7f, absorb_red_blue);
  Material fuzzy_mirror = Material::Mirror(0.2);

  // Background lighting
  scene.SetAmbientBackground(Spectrum::White * 0.01);
  scene.SetSkylight(Spectrum::Cyan * 0.01, -up);

  // Ground plane
  Vec3 O1 = kZero3;
  scene.Add(Plane(O1, O1 + right, O1 + forward, lambertian));

  // Big sphere
  float r2 = 0.9f;
  Vec3 O2 = O1 + 2.5 * forward + r2 * up;
  scene.Add(Sphere(O2, r2, fuzzy_mirror));

  // Small green sphere on right
  float r3 = 0.4f;
  Vec3 O3 = O2 - 1.5 * forward + 0.5 * right - r2 * up + r3 * up;
  scene.Add(Sphere(O3, r3, clean_green_metal));

  // Medium blue sphere on left
  float r4 = 0.5f;
  Vec3 O4 = O2 - 1.3 * forward - 0.5 * right - r2 * up + r4 * up;
  scene.Add(Sphere(O4, r4, clean_blue_metal));

  // Camera
  Vec3 Ocam = O1 + 0.6 * up;
  Camera cam(Ocam, (O2 - r2 * up - Ocam).NormalizedNonzero(), right, width,
             height);
  scene.Add(cam);

  return scene;
}

// A static scene inside a box with reflective sides, some of which are tinted.
// The box contains a metallic sphere and a small spherical light source.
// Scene arguments: none
Scene mirror_box(int width, int height, arglist args = {}) {
  Scene scene;

  // Camera
  Vec3 cam_pos = kZero3;
  Camera cam(cam_pos, forward, right, width, height);
  scene.Add(cam);

  // Materials
  Material white_light = Material::Light();
  Material dim_white_light = Material::Light(Spectrum::White * 0.05);
  Material metal = Material::Metal(0.8f);
  Material red_metal = Material::Metal(0.8f, Spectrum(Gaussian(1.0f, 500, 80)));

  // Auxiliary screen & direction variables
  float aspect = static_cast<float>(width) / height;
  Vec3 screen_top = (right ^ forward) / aspect;

  // Big metal sphere
  Vec3 sphere_vel(0, 0, 0);
  Vec3 sphere_O = 2 * forward;
  float sphere_rad = 1.0f;
  Sphere sphere(sphere_O, sphere_rad, sphere_vel, 0, metal);
  scene.Add(sphere);

  // Small spherical light
  Vec3 light_sphere_O = cam_pos + 0.5 * right + 1.5f * screen_top;
  float light_sphere_rad = 0.6f;
  Sphere light_sphere(light_sphere_O, light_sphere_rad, white_light);
  scene.Add(light_sphere);

  // Surrounding box - tint right side red, make top emit dim light
  Vec3 box_O = sphere_O;
  float box_side = 3.0f;
  Box surrounding_box(box_O, right * box_side, forward * box_side,
                      up * box_side);
  surrounding_box.back.SetMaterial(metal);
  surrounding_box.down.SetMaterial(metal);
  surrounding_box.right.SetMaterial(red_metal);
  surrounding_box.up.SetMaterial(dim_white_light);
  surrounding_box.left.SetMaterial(metal);
  surrounding_box.front.SetMaterial(metal);
  scene.Add(surrounding_box);

  return scene;
}

// A square box moving right. Scene arguments allow customisation of speed and
// material properties of the box. Arguments from the fifth onwards are
// interpreted as data of absorption modes of the box.
// Scene arguments: speed albedo diffuse specular fuzz amplitude_1 mean_1
//  sigma_1 amplitude_2 mean_2 sigma_2 amplitude_3 ...
Scene moving_box(int width, int height, arglist args = {}) {
  // Standard box properties - customisable by scene arguments.
  float speed = 0.25f;
  float albedo = 0.8f;
  float diffuse = 0.0f;
  float specular = 1.0f;
  float fuzz = 0.2;
  Spectrum absorption = Spectrum::Black;

  // Read in scene arguments provided
  if (args.size() > 0) {
    speed = std::stof(args[0]);
  }
  if (args.size() > 1) {
    albedo = std::stof(args[1]);
  }
  if (args.size() > 2) {
    diffuse = std::stof(args[2]);
  }
  if (args.size() > 3) {
    specular = std::stof(args[3]);
  }
  if (args.size() > 4) {
    fuzz = std::stof(args[4]);
  }
  if (args.size() > 5) {
    // Read in absorption modes
    int i = 5;
    while (i + 2 < args.size()) {
      float amplitude = std::stof(args[i]);
      float lambda = std::stof(args[i + 1]);
      float sigma = std::stof(args[i + 2]);
      i += 3;
    }
    if (i < args.size()) {
      std::cout << "Ignoring remaining scene arguments, "
                << "starting with" << args[i] << std::endl;
    }
  }
  std::cout << "Setting up moving box with "
            << "speed " << speed << " "
            << "albedo " << albedo << " "
            << "diffuse " << diffuse << " "
            << "specular " << specular << " "
            << "fuzz " << fuzz << std::endl;
  if (absorption.size() < 1) {
    std::cout << "and no absorption." << std::endl;
  } else {
    std::cout << "and absorption modes " << std::endl;
    for (auto mode : absorption.getModes()) {
      std::cout << "(" << mode.amplitude << ", " << mode.mean << ", "
                << mode.sig << ")" << std::endl;
    }
  }

  // Finished reading in scene arguments. Begin setting up scene.
  Scene scene;
  scene.SetAmbientBackground(Spectrum::White * 0.01);

  // Materials
  Material clean_metal = Material::Metal(0.8f);
  Material box_material = Material(albedo, diffuse, specular, fuzz, absorption);

  // Box
  Vec3 box_vel = speed * right;
  Vec3 box_O = 3 * forward;
  float box_side = 0.5f;
  Box box(box_O, right * box_side, forward * box_side, up * box_side, box_vel,
          0, box_material);
  scene.Add(box);

  // Flat plane
  Vec3 plane_O = box_O + down * 3 * box_side / 2.0f;
  Vec3 plane_a = right;
  Vec3 plane_b = forward;
  Plane plane(plane_O, plane_a, plane_b, clean_metal);
  scene.Add(plane);

  // Camera
  Vec3 cam_pos = 3 * box_side * up;
  Camera cam(cam_pos, (box_O - cam_pos).NormalizedNonzero(), right, width,
             height);
  scene.Add(cam);

  return scene;
}

// A row of boxes moving up/down at varying speeds above a mirror.
// Fast-moving boxes are very distorted. The images in the mirror are
// time-delayed because of longer light travel time.
// The number of boxes and min/max speeds are customizable. Providing the flag
// "bg" in the 4th position activates a vertical background plane.
// Scene arguments: nr_boxes min_speed max_speed bg
Scene boxes_in_mirror(int width, int height, arglist args = {}) {
  // Standard scene properties
  int nr_boxes = 5;
  float min_speed = -0.5f;
  float max_speed = 0.5f;
  bool bg_activated = false;

  // Read in scene arguments provided
  int nargs = args.size();
  if (nargs > 0) {
    nr_boxes = std::stoi(args[0]);
  }
  if (nargs > 1) {
    min_speed = std::stof(args[1]);
  }
  if (nargs > 2) {
    max_speed = std::stof(args[2]);
  }
  if (nargs > 3 && args[3] == "bg") {
    bg_activated = true;
  }

  // Finished reading scene args. Set up scene.
  Scene scene;
  scene.SetAmbientBackground(Spectrum::White);

  // Camera data
  Vec3 cam_pos = kZero3;
  Vec3 cam_to_screen = forward + 0.2 * down;
  float angle = 60 * kPi / 180;  // vertical screen angle
  float screen_height = 2 * cam_to_screen.Norm() * tanf(angle / 2);
  float screen_width = width * screen_height / height;

  // Box data
  float d_speed = (max_speed - min_speed) / (nr_boxes - 1);
  float side = 1.0f;            // sidelength of boxes
  float spacing = 0.5f * side;  // spacing between boxes
  Material box_material = Material(0.5);
  float distance_to_cam =
      nr_boxes * (side + spacing) * width / (2.f * tanf(angle / 2.f) * height);
  Vec3 delta_centre = (side + spacing) * right;
  Vec3 leftmost_centre =
      cam_pos + forward * distance_to_cam - (nr_boxes - 1) / 2.f * delta_centre;
  Vec3 c2r = side / 2.f * right;
  Vec3 c2b = side / 2.f * forward;
  Vec3 c2u = side / 2.f * up;

  // Set up boxes. Arrange them such that at t=0, their centers are all aligned
  // in the middle row of the screen.
  for (int box_index = 0; box_index < nr_boxes; box_index++) {
    float speed = min_speed + box_index * d_speed;
    Vec3 centre = leftmost_centre + delta_centre * box_index;
    // travel time for a lightray to get from the box centre to the camera
    float time_offset = (cam_pos - centre).Norm();
    Box box(centre, c2r, c2b, c2u, speed * up, -time_offset, box_material);
    scene.Add(box);
  }

  // Add camera
  Camera cam(cam_pos, cam_to_screen, screen_width / 2.f * right, width, height);
  scene.Add(cam);

  // Add ground mirror plane
  Material mirror_material = Material::Metal(0.8f);
  float mirror_distance_to_cam = 2 * side;
  Plane mirror(cam_pos + mirror_distance_to_cam * down, right, forward,
               mirror_material);
  scene.Add(mirror);

  // Add vertical background plane (if activated)
  if (bg_activated) {
    Material bg_material = Material(0.01f, 0.8f, 0.2f);
    float bg_distance_to_cam = distance_to_cam + side * 0.7f;
    Plane bg_plane(
        cam_pos + bg_distance_to_cam * cam_to_screen.NormalizedNonzero(), right,
        up, bg_material);
    scene.Add(bg_plane);
  }

  return scene;
}

// Three identical spheres emitting green light are placed in front of a mirror.
// The top/bottom sphere is moving towards/away from the camera and away/towards
// the mirror. This scene shows the Doppler shift and time delay due to light
// travel times. Specifying scene arguments allows changing the number of
// spheres as well as the magnitude of the speeds.
// Scene arguments: nr_spheres max_speed
Scene glowing_spheres(int width, int height, arglist args = {}) {
  // Standard scene arguments
  float max_speed = 0.25f;
  int nr_spheres = 3;
  if (args.size() > 0) {
    max_speed = std::stof(args[0]);
  }
  if (args.size() > 1) {
    nr_spheres = std::stoi(args[1]);
  }

  // Finished reading in scene arguments. Set up scene.
  Scene scene;
  scene.SetAmbientBackground(Spectrum::White * 0.01);

  // Camera
  Vec3 cam_pos = kZero3;
  Camera cam(cam_pos, forward, right, width, height);
  scene.Add(cam);

  // Materials
  Material tinted_mirror = Material::Metal(0.6);
  Material white_light =
      Material::Metal().SetProperty("emission", Spectrum::White);
  Spectrum green = Spectrum(Gaussian(0.5f, 530, 20));
  Material green_glowing = Material::Metal(0.0f).SetProperty("emission", green);

  // Sphere and mirror data
  float sphere_distance = 3.5f;
  float offset_angle = 15 * kPi / 180;
  float mirror_distance = sphere_distance * 1.9f;
  float mirror_angle = offset_angle * 2;

  // Add spheres
  float sphere_rad = 0.7f;
  float spacing = sphere_rad;
  Vec3 delta_centre = (2.0 * sphere_rad + spacing) * down;
  Vec3 top_sphere_centre =
      (forward * cosf(offset_angle) + right * sinf(offset_angle)) *
          sphere_distance -
      delta_centre * (nr_spheres - 1) / 2.0f;
  float min_speed = -max_speed;
  float delta_speed = 2 * max_speed / (nr_spheres - 1);
  for (int i = 0; i < nr_spheres; i++) {
    Vec3 sphere_vel = (min_speed + i * delta_speed) * forward;
    Vec3 sphere_O = top_sphere_centre + delta_centre * i;
    float time_correction = (sphere_O - cam_pos).Norm();
    Sphere sphere(sphere_O, sphere_rad, sphere_vel, -time_correction,
                  green_glowing);
    scene.Add(sphere);
  }

  // Add tinted mirror
  Vec3 mirror_O = cam_pos + forward * mirror_distance;
  Vec3 mirror_A = (right * cosf(mirror_angle) + forward * sinf(mirror_angle)) *
                  4 * sphere_rad;
  Vec3 mirror_B = up * nr_spheres * (2 * sphere_rad + spacing) / 2.f;
  Parallelogram mirror(mirror_O - mirror_A / 2.0f, mirror_A, mirror_B,
                       tinted_mirror);
  scene.Add(mirror);

  return scene;
}

// Seven boxes arranged vertically are moving sideways with differing speeds.
// The box material has a single absorption mode in the short wavelength range
// (amplitude: 1, mean: 400nm, std: 50nm). The scene shows the effect of the
// Doppler shift on such an absorption mode. The absorption properties
// may be changed by specifying absorption mode data as scene arguments.
// Scene arguments: amplitude_1 mean_1 sigma_1 amplitude_2 ...
Scene box_array(int width, int height, arglist args = {}) {
  Spectrum box_absorption(Gaussian(1.0f, 400, 50));  // standard absorption

  // Read in scene arguments
  int n_args = args.size();
  if (n_args > 0) {
    Spectrum new_absorb = Spectrum::Black;
    int i = 0;
    while (i + 2 < n_args) {
      float amp = std::stof(args[i]);
      float mean = std::stof(args[i + 1]);
      float sig = std::stof(args[i + 2]);
      new_absorb += Gaussian(amp, mean, sig);
      i += 3;
    }
    if (i < n_args) {
      std::cout << "Ignoring remaining arguments, "
                << "starting at " << args[i] << std::endl;
    }
    box_absorption = new_absorb;
  }

  // Create box material with chosen absorption properties.
  Material box_material =
      Material::Metal(0.8f);  // let this be changed by arguments?
  box_material.SetProperty("absorption", box_absorption);

  // Finished dealing with scene arguments. Set up scene.
  Scene scene;
  scene.SetAmbientBackground(Spectrum::White);

  // Box data
  int nr_boxes = 7;
  float min_speed = -0.7f;
  float max_speed = 0.7f;
  float d_speed = (max_speed - min_speed) / (nr_boxes - 1);
  float side = 1.0f;            // sidelength of boxes
  float spacing = 1.5f * side;  // spacing between boxes

  // Camera
  Vec3 cam_pos = kZero3;
  Vec3 cam_to_screen = forward;
  float angle = 60 * kPi / 180;  // vertical screen angle
  float screen_height = 2 * cam_to_screen.Norm() * tanf(angle / 2);
  float screen_width = width * screen_height / height;
  Camera cam(cam_pos, cam_to_screen, screen_width / 2.f * right, width, height);
  scene.Add(cam);

  // Compute box data
  float distance_to_cam =
      nr_boxes * (side + spacing) / (2.f * tanf(angle / 2.f));
  Vec3 lowest_centre = (side + spacing) * (nr_boxes - 1) / 2.f * down +
                       distance_to_cam * cam_to_screen.NormalizedNonzero();
  Vec3 delta_centre = (side + spacing) * up;
  Vec3 c2r = side / 2.f * right;
  Vec3 c2b = side / 2.f * forward;
  Vec3 c2u = side / 2.f * up;

  // Add boxes. Arrange them such that at t=0, their centers are all aligned
  // in the central column of the screen.
  for (int box_index = 0; box_index < nr_boxes; box_index++) {
    float speed = min_speed + box_index * d_speed;
    Vec3 center = lowest_centre + delta_centre * box_index;
    // travel time for a lightray to get from the box center to the camera
    float time_offset = (cam_pos - center).Norm();
    Box box(center, c2r, c2b, c2u, speed * right, -time_offset, box_material);
    scene.Add(box);
  }

  // Add vertical background plane behind the boxes.
  float plane_distance_to_cam = 2 * distance_to_cam;
  Plane plane(
      cam_pos + plane_distance_to_cam * cam_to_screen.NormalizedNonzero(),
      right, up, Material::Metal(0.1, Spectrum::Black, 0.8f));
  scene.Add(plane);

  return scene;
}

// Two identical balls travel downwards at the same speed. The view towards
// the left ball is obstructed by a glass box (top) and a diamond box (bottom).
// This scene shows the time delay resulting from slower light travel time in
// dielectric media, as well as the Doppler effect, visible especially for the
// reflections of the balls in the sides of the boxes. Tjhe speed of the boxes
// can be changed by specifying it as a scene argument.
// Scene arguments: ball_speed
Scene dielectric_delays(int width, int height, arglist args = {}) {
  float ball_speed = 0.25f;
  if (args.size() > 0) {
    ball_speed = std::stof(args[0]);
  }

  // Finished reading in scene argument. Set up scene.
  Scene scene;
  scene.SetSkylight(Spectrum::White * 0.05,
                    (down + forward + left).NormalizedNonzero());
  scene.SetAmbientBackground(Spectrum::White * 0.1);

  // Materials
  Material glass = Material::Glass;
  Material diamond = Material::Diamond;

  // Camera
  float horizontal_fov = 60 * kPi / 180;
  float half_fov = horizontal_fov / 2;
  float quarter_fov = half_fov / 2;
  Camera cam(kZero3, forward, right * tanf(half_fov), width, height);
  scene.Add(cam);

  // Add boxes
  float box_height = 1.0f;
  float box_width = 1.0f;
  float box_depth = 1.0f;
  float spacing = 0.5f;
  float mid_distance = 3.0f + box_depth / 2;
  Vec3 mid_left = mid_distance * (forward + tanf(quarter_fov) * left);
  Vec3 g_center = mid_left + up * (box_height + spacing) / 2;
  Vec3 d_center = mid_left + down * (box_height + spacing) / 2;
  Vec3 c2r = right * box_width / 2;
  Vec3 c2b = forward * box_depth / 2;
  Vec3 c2u = up * box_height / 2;
  Box glass_box(g_center, c2r, c2b, c2u, glass);
  Box diamond_box(d_center, c2r, c2b, c2u, diamond);
  scene.Add(glass_box);
  scene.Add(diamond_box);

  // Add travelling balls
  float far_distance = mid_distance + box_depth;
  Vec3 far_left = far_distance * (forward + tanf(quarter_fov) * left);
  Vec3 far_right = far_distance * (forward + tanf(quarter_fov) * right);
  Vec3 ball_vel = ball_speed * down;
  float ball_rad = 0.2 * box_height;
  Material blue_emitter = Material::Light(Spectrum(Gaussian(1.0f, 400, 50)));
  Sphere left_ball(far_left, ball_rad, ball_vel, 0.0f, blue_emitter);
  Sphere right_ball(far_right, ball_rad, ball_vel, 0.0f, blue_emitter);
  scene.Add(left_ball);
  scene.Add(right_ball);

  return scene;
}

// Headlight absorption effect inside cube. Three identical, diffusely
// reflecting, white balls are placed inside a box. The top ball moves right,
// the bottom ball to the light, the middle ball is stationary. The left/right
// panel of the box emits red/green light. This scene shows the headlight effect
// for absorption. The speed of the moving boxes can be changed via a scene
// argument.
// Scene arguments: speed
Scene headlight_absorption(int width, int height, arglist args = {}) {
  float speed = 0.4f;
  if (args.size() > 0) {
    speed = std::stof(args[0]);
  }

  // Finished reading in scene argument. Set up scene.
  Scene scene;

  // Camera
  Vec3 cam_pos(kZero3);
  float horz_fov = 100 * kPi / 180;
  float half_fov = horz_fov / 2;
  Camera cam(cam_pos, forward, right * tanf(half_fov), width, height);
  scene.Add(cam);

  // Add surrounding box, replace side panels material with lights.
  // The remaining four panels are set to a fuzzy metal.
  float box_side = 5.0f;
  float half_side = box_side / 2;
  Vec3 box_O = cam_pos + forward * box_side * 0.4f;
  Box box(box_O, right * half_side, forward * half_side, up * half_side,
          Material::Metal(0.5f, Spectrum::Black, 0.5f));
  box.right.SetMaterial(Material::Light(Spectrum(Gaussian(1.0f, 625, 25))));
  box.left.SetMaterial(Material::Light(Spectrum(Gaussian(1.0f, 525, 25))));
  scene.Add(box);

  // Add three balls
  float ball_rad = box_side / 10.f;
  Material lambertian(0.95f);
  Vec3 mid_ball_pos = box_O + forward * box_side / 4;
  Vec3 top_ball_pos = mid_ball_pos + up * 3 * ball_rad;
  Vec3 bot_ball_pos = mid_ball_pos + down * 3 * ball_rad;
  Sphere mid_ball(mid_ball_pos, ball_rad, lambertian);
  float time_delay = (top_ball_pos - cam_pos).Norm() + 2 * ball_rad;
  Sphere top_ball(top_ball_pos, ball_rad, speed * right, -time_delay,
                  lambertian);
  Sphere bot_ball(bot_ball_pos, ball_rad, speed * left, -time_delay,
                  lambertian);
  scene.Add(top_ball);
  scene.Add(mid_ball);
  scene.Add(bot_ball);

  return scene;
}

// An empty custom scene.
Scene custom_scene(int width, int height, arglist args = {}) {
  Scene scene;

  return scene;
}

const std::map<std::string, Scene (*)(int, int, arglist)> scene_map = {
    {"glowing_spheres", &glowing_spheres},
    {"moving_box", &moving_box},
    {"boxes_in_mirror", &boxes_in_mirror},
    {"metal_balls", &metal_balls},
    {"box_array", &box_array},
    {"mirror_box", &mirror_box},
    {"dielectric_delays", &dielectric_delays},
    {"headlight_absorption", &headlight_absorption},
    {"custom_scene", &custom_scene}};

int main(int argc, char *argv[]) {
  // Start with standard arguments, then parse command line args.
  cmd_line_args args;
  if (argc > 1) {
    bool parsed = args.parse_arguments(argc, argv);
    if (!parsed) {
      std::cout << "Error parsing cmd line arguments. Terminating."
                << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Successfully parsed cmd line arguments." << std::endl;
  }

  // Do tests if specified
  if (args.test) {
    Tester test;
    bool success = false;
    if (args.scene_name == "lorentz") {
      success =
          test.TestLineLorentzTransforms() && test.TestColorLorentzTransforms();
    } else if (args.scene_name == "line_lorentz") {
      success = test.TestLineLorentzTransforms();
    } else if (args.scene_name == "color_lorentz") {
      success = test.TestColorLorentzTransforms();
    } else if (args.scene_name == "color") {
      success = test.TestColors();
    } else if (args.scene_name == "scenes") {
      success = test.TestScenes();
    } else {
      // default: run all tests
      success = test.RunAllTests();
    }

    if (success) return EXIT_SUCCESS;
    return EXIT_FAILURE;
  }

  // Render the specified scene:
  auto scene_iter = scene_map.find(args.scene_name);
  if (scene_iter == scene_map.end()) {
    std::cout << "Specified scene not found. Terminating." << std::endl;
    return EXIT_FAILURE;
  }  // else
  auto scene_fn = scene_iter->second;
  Scene scene = scene_fn(args.width, args.height, args.scene_args);

  // Raytracer
  Tracer tracey(scene);
  if (!args.film) {
    tracey.RenderImage(args.rays_per_pixel, args.depth, args.camera_index,
                       args.iterations_per_update, args.filename, args.time,
                       args.rescale_factor);
  } else {
    tracey.RenderFilm(args.rays_per_pixel, args.depth, args.camera_index,
                      args.iterations_per_update, args.filename, args.time,
                      args.end_time, args.d_time, args.start_frame,
                      args.end_frame, args.preview_only, args.rescale_factor);
  }

  return EXIT_SUCCESS;
}
