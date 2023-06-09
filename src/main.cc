// Copyright 2023 Christoph Weis
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>

#include "include/test.h"
#include "include/tracer.h"

typedef std::vector<std::string> arglist;

struct cmd_line_args {
  // Scene and output args
  std::string scene_name = "box_array";
  std::string filename = "output";
  bool film = false;
  bool preview_only = false;
  bool test = false;  // run tests instead of render

  // Standard rendering arguments
  int width = 500;
  int height = 500;
  int rays_per_pixel = 100;
  int depth = 8;
  int iterations_per_update = 1;
  int camera_index = 0;

  // Time arguments
  float time = -1.0f;     // = start_time for film
  float end_time = 3.0f;  // ignored for image
  float d_time = 0.1f;    // ignored for image

  // Parallel processing helps
  int start_frame = 0;        // first frame to render
  int end_frame = -1;         // code for: do all frames
  float rescale_factor = -1;  // code for: estimate rescale factor

  // Scene arguments
  arglist scene_args = {};

  // Now split up arguments by type and give names and aliases.
  std::map<std::string, int &> int_args = {
      {"--width", width},
      {"-w", width},
      {"--height", height},
      {"-h", height},
      {"--rays_per_pixel", rays_per_pixel},
      {"-rpp", rays_per_pixel},
      {"--depth", depth},
      {"-d", depth},
      {"--camera_index", camera_index},
      {"-ci", camera_index},
      {"--iterations_per_update", iterations_per_update},
      {"-ipu", iterations_per_update},
      {"--start_frame", start_frame},
      {"-sf", start_frame},
      {"--end_frame", end_frame},
      {"-ef", end_frame}};

  std::map<std::string, float &> float_args = {
      {"--start_time", time},
      {"-st", time},
      {"--time", time},
      {"-t", time},
      {"--end_time", end_time},
      {"-et", end_time},
      {"--d_time", d_time},
      {"-dt", d_time},
      {"--rescale_factor", rescale_factor},
      {"-rf", rescale_factor}};

  std::map<std::string, std::string &> string_args = {
      {"--scene_name", scene_name},
      {"-sn", scene_name},
      {"--filename", filename},
      {"-fn", filename}};

  std::map<std::string, bool &> bool_args = {{"--test", test},
                                             {"-t", test},
                                             {"--film", film},
                                             {"-f", film},
                                             {"--preview_only", preview_only},
                                             {"-po", preview_only}};

  std::vector<std::string> scene_args_triggers = {"--scene_args", "--args",
                                                  "-sa"};
};

// Parses command line arguments
bool parse_arguments(int argc, char *argv[], cmd_line_args &args) {
  // skip first argument: it's the filename
  for (int i = 1; i < argc; i++) {
    // check if it's an integer argument
    {
      auto int_iter = args.int_args.find(argv[i]);
      if (int_iter != args.int_args.end()) {
        if (i + 1 >= argc) {
          std::cout << "Parameter keyword '" << argv[i]
                    << "' specified without value." << std::endl;
          return false;
        }
        int_iter->second = std::stoi(argv[i + 1]);
        i++;
        continue;
      }
    }

    // check if it's a float argument
    {
      auto float_iter = args.float_args.find(argv[i]);
      if (float_iter != args.float_args.end()) {
        if (i + 1 >= argc) {
          std::cout << "Parameter keyword '" << argv[i]
                    << "' specified without value." << std::endl;
          return false;
        }
        float_iter->second = std::stof(argv[i + 1]);
        i++;
        continue;
      }
    }

    // check if it's a string argument
    {
      auto string_iter = args.string_args.find(argv[i]);
      if (string_iter != args.string_args.end()) {
        if (i + 1 >= argc) {
          std::cout << "Parameter keyword '" << argv[i]
                    << "' specified without value." << std::endl;
          return false;
        }
        string_iter->second = argv[i + 1];
        i++;
        continue;
      }
    }

    // check if it's a bool argument
    {
      auto bool_iter = args.bool_args.find(argv[i]);
      if (bool_iter != args.bool_args.end()) {
        bool_iter->second = true;
        continue;
      }
    }

    // check if it's the start of the scene arguments
    {
      auto sa_trigger_iter =
          std::find(args.scene_args_triggers.cbegin(),
                    args.scene_args_triggers.cend(), argv[i]);
      if (sa_trigger_iter < args.scene_args_triggers.cend()) {
        // add all following arguments into the scene_args vector
        for (i = i + 1; i < argc; i++) {
          args.scene_args.push_back(argv[i]);
        }
        break;
      }
    }

    // it didn't match any of the above
    std::cout << "Ignoring argument '" << argv[i]
              << "', as it doesn't match any parameters." << std::endl;
  }

  return true;
}

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

Scene metal_balls(int width, int height, arglist args = {}) {
  Scene scene;

  Vec3 forward = y;
  Vec3 right = x;
  Vec3 up = right ^ forward;

  Material lambertian = Material(0.9f, 0.7f);
  Spectrum absorb_red_green(
      {GaussianData(0.5f, 550, 50), GaussianData(0.5f, 600, 50)});
  Spectrum absorb_red_blue(
      {GaussianData(0.5f, 450, 50), GaussianData(0.5f, 600, 25)});
  Material clean_blue_metal = Metal(0.7f, absorb_red_green);
  Material clean_green_metal = Metal(0.7f, absorb_red_blue);
  Material fuzzy_mirror = Mirror(0.2);

  scene.SetAmbientBackground(kWhite * 0.01);
  scene.SetSkylight(kDimCyan * 0.1, -up);

  Vec3 O1 = kZero3;
  scene.Add(Plane(O1, O1 + right, O1 + forward, lambertian));

  float r2 = 0.9f;
  Vec3 O2 = O1 + 2.5 * forward + r2 * up;
  scene.Add(Sphere(O2, r2, fuzzy_mirror));

  float r3 = 0.4f;
  Vec3 O3 = O2 - 1.5 * forward + 0.5 * right - r2 * up + r3 * up;
  scene.Add(Sphere(O3, r3, clean_green_metal));

  float r4 = 0.5f;
  Vec3 O4 = O2 - 1.3 * forward - 0.5 * right - r2 * up + r4 * up;
  scene.Add(Sphere(O4, r4, clean_blue_metal));

  Vec3 Ocam = O1 + 0.6 * up;
  Camera cam(Ocam, (O2 - r2 * up - Ocam).NormalizedNonzero(), right, width,
             height);
  scene.AddCamera(cam);

  return scene;
}

Scene moving_box(int width, int height, arglist args = {}) {
  /*
    MOVING BOX
  */

  // Scene args - standard values
  //  speed & material properties
  //  of box
  float speed = 0.25f;
  float albedo = 0.8f;
  float diffuse = 0.0f;
  float specular = 1.0f;
  float fuzz = 0.2;
  Spectrum absorption = kBlack;

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
  if (absorption.getModes().size() <= 5) {
    std::cout << "and no absorption." << std::endl;
  } else {
    std::cout << "and absorption modes " << std::endl;
    for (auto mode : absorption.getModes()) {
      std::cout << "(" << mode.amplitude << ", " << mode.mean << ", "
                << mode.sig << ")" << std::endl;
    }
  }

  Scene scene;

  // auxiliary
  Vec3 forward = z;
  Vec3 right = x;
  Vec3 up = Cross(right, forward);
  Vec3 backward = -forward;
  Vec3 left = -right;
  Vec3 down = -up;

  scene.SetAmbientBackground(kWhite * 0.01);

  // Materials
  Material clean_metal = Metal(0.8f, kBlack);
  Material box_material = Material(albedo, diffuse, specular, fuzz, absorption);

  // Add objects
  // Vec3 sphere_vel = 0.9 * right;
  Vec3 box_vel = speed * right;
  Vec3 box_O = 3 * forward;
  float box_side = 0.5f;

  Box box(box_O, right * box_side, forward * box_side, up * box_side, box_vel,
          0, box_material);
  scene.Add(box);

  Vec3 plane_O = box_O + down * 3 * box_side / 2.0f;
  Vec3 plane_a = right;
  Vec3 plane_b = forward;
  Plane plane(plane_O, plane_a, plane_b, clean_metal);
  scene.Add(plane);

  Vec3 cam_pos = 3 * box_side * up;
  Camera cam(cam_pos, (box_O - cam_pos).NormalizedNonzero(), right, width,
             height);
  scene.AddCamera(cam);

  return scene;
}

Scene mirror_box(int width, int height, arglist args = {}) {
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

  Material white_light = Light();
  Material dim_white_light = Light(kWhite * 0.05);
  Material metal = Metal(0.8f);
  Material red_metal = Metal(0.8f, Spectrum(GaussianData(1.0f, 500, 80)));
  Material glass = kGlass;

  // Auxiliary screen & direction variables
  float aspect = static_cast<float>(width) / height;
  Vec3 screen_top = (right ^ forward) / aspect;

  // Add objects
  Vec3 sphere_vel(0, 0, 0);
  Vec3 sphere_O = 2 * forward;
  float sphere_rad = 1.0f;
  Sphere sphere(sphere_O, sphere_rad, sphere_vel, 0, metal);
  scene.Add(sphere);

  Vec3 light_sphere_O = cam_pos + 0.5 * right + 1.5f * screen_top;
  float light_sphere_rad = 0.6f;
  Sphere light_sphere(light_sphere_O, light_sphere_rad, white_light);
  scene.Add(light_sphere);

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

  std::cout << "Scene created." << std::endl << std::flush;

  return scene;
}

Scene boxes_in_mirror(int width, int height, arglist args = {}) {
  int nr_boxes = 5;
  float min_speed = -0.5f;
  float max_speed = 0.5f;
  bool bg_activated = false;

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

  float d_speed = (max_speed - min_speed) / (nr_boxes - 1);

  float side = 1.0f;            // sidelength of boxes
  float spacing = 0.5f * side;  // spacing between boxes

  Material box_material = Material(0.5);

  Vec3 up(0, 1, 0);
  Vec3 right(1, 0, 0);
  Vec3 forward = Cross(up, right);
  Vec3 down = -up;
  Vec3 left = -right;
  Vec3 backward = -forward;

  Scene scene;
  scene.SetAmbientBackground(kWhite);

  Vec3 cam_pos = kZero3;
  Vec3 cam_to_screen = forward + 0.2 * down;
  float angle = 60 * kPi / 180;  // vertical screen angle
  float screen_height = 2 * cam_to_screen.Norm() * tanf(angle / 2);
  float screen_width = width * screen_height / height;

  float distance_to_cam =
      nr_boxes * (side + spacing) * width / (2.f * tanf(angle / 2.f) * height);
  Vec3 delta_centre = (side + spacing) * right;
  Vec3 leftmost_centre =
      cam_pos + forward * distance_to_cam - (nr_boxes - 1) / 2.f * delta_centre;
  Vec3 c2r = side / 2.f * right;
  Vec3 c2b = side / 2.f * forward;
  Vec3 c2u = side / 2.f * up;

  // arrange boxes such that at t=0, they are all aligned in the middle row
  // of the screen
  for (int box_index = 0; box_index < nr_boxes; box_index++) {
    float speed = min_speed + box_index * d_speed;
    Vec3 centre = leftmost_centre + delta_centre * box_index;
    // travel time for a lightray to get from the box centre to the camera
    float time_offset = (cam_pos - centre).Norm();
    Box box(centre, c2r, c2b, c2u, speed * up, -time_offset, box_material);
    scene.Add(box);
  }

  Camera cam(cam_pos, cam_to_screen, screen_width / 2.f * right, width, height);
  scene.AddCamera(cam);

  Material mirror_material = Metal(0.8f);
  float mirror_distance_to_cam = 2 * side;
  Plane mirror(cam_pos + mirror_distance_to_cam * down, right, forward,
               mirror_material);
  scene.Add(mirror);

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

Scene glowing_spheres(int width, int height, arglist args = {}) {
  /*
    MOVING SPHERE
  */
  float max_speed = 0.25f;  // standard value
  int nr_spheres = 3;
  if (args.size() > 0) {
    max_speed = std::stof(args[0]);
  }
  if (args.size() > 1) {
    nr_spheres = std::stoi(args[1]);
  }

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
  Material mirror = Metal(0.6);
  Material white_light = Metal().SetProperty("emission", kWhite);
  Spectrum green = Spectrum(GaussianData(0.5f, 530, 20));
  Material green_glowing = Metal(0.0f).SetProperty("emission", green);

  float sphere_distance = 3.5f;
  float offset_angle = 15 * kPi / 180;
  float mirror_distance = sphere_distance * 1.9f;
  float mirror_angle = offset_angle * 2;

  // Add objects
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

  Vec3 mirror_O = cam_pos + forward * mirror_distance;
  Vec3 mirror_A = (right * cosf(mirror_angle) + forward * sinf(mirror_angle)) *
                  4 * sphere_rad;
  Vec3 mirror_B = up * nr_spheres * (2 * sphere_rad + spacing) / 2.f;
  Parallelogram mirror_plane(mirror_O - mirror_A / 2.0f, mirror_A, mirror_B,
                             mirror);
  scene.Add(mirror_plane);

  scene.SetAmbientBackground(kWhite * 0.01);

  return scene;
}

// Stack boxes underneath each other, with increasing speed to the left, going
// from bottom to top.
Scene box_array(int width, int height, arglist args = {}) {
  int nr_boxes = 7;
  float min_speed = -0.7f;
  float max_speed = 0.7f;
  // float d_speed = 0.1f;
  float d_speed = (max_speed - min_speed) / (nr_boxes - 1);

  float side = 1.0f;            // sidelength of boxes
  float spacing = 1.5f * side;  // spacing between boxes

  Spectrum absorb_short(GaussianData(1.0f, 400, 50));

  int n_args = args.size();
  if (n_args > 0) {
    Spectrum new_absorb = kBlack;
    int i = 0;
    while (i + 2 < n_args) {
      float amp = std::stof(args[i]);
      float mean = std::stof(args[i + 1]);
      float sig = std::stof(args[i + 2]);
      new_absorb += GaussianData(amp, mean, sig);
      i += 3;
    }

    if (i < n_args) {
      std::cout << "Ignoring remaining arguments, "
                << "starting at " << args[i] << std::endl;
    }
    absorb_short = new_absorb;
  }

  Material box_material = Metal(0.8f);  // let this be changed by arguments?
  box_material.SetProperty("absorption", absorb_short);

  Vec3 up(0, 1, 0);
  Vec3 right(1, 0, 0);
  Vec3 forward = Cross(up, right);
  Vec3 down = -up;
  Vec3 left = -right;
  Vec3 backward = -forward;

  Scene scene;
  scene.SetAmbientBackground(kWhite);
  // scene.SetSkylight(kWhite, forward + down);

  Vec3 cam_pos = kZero3;
  Vec3 cam_to_screen = forward;
  float angle = 60 * kPi / 180;  // vertical screen angle
  float screen_height = 2 * cam_to_screen.Norm() * tanf(angle / 2);
  float screen_width = width * screen_height / height;
  Camera cam(cam_pos, cam_to_screen, screen_width / 2.f * right, width, height);
  scene.AddCamera(cam);

  float distance_to_cam =
      nr_boxes * (side + spacing) / (2.f * tanf(angle / 2.f));
  Vec3 lowest_centre = (side + spacing) * (nr_boxes - 1) / 2.f * down +
                       distance_to_cam * cam_to_screen.NormalizedNonzero();
  Vec3 delta_centre = (side + spacing) * up;
  Vec3 c2r = side / 2.f * right;
  Vec3 c2b = side / 2.f * forward;
  Vec3 c2u = side / 2.f * up;

  // arrange boxes such that at t=0, they are all aligned in the middle column
  // of the screen
  for (int box_index = 0; box_index < nr_boxes; box_index++) {
    float speed = min_speed + box_index * d_speed;
    Vec3 centre = lowest_centre + delta_centre * box_index;
    // travel time for a lightray to get from the box centre to the camera
    float time_offset = (cam_pos - centre).Norm();
    Box box(centre, c2r, c2b, c2u, speed * right, -time_offset, box_material);
    scene.Add(box);
  }

  float plane_distance_to_cam = 2 * distance_to_cam;
  Plane plane(
      cam_pos + plane_distance_to_cam * cam_to_screen.NormalizedNonzero(),
      right, up, Metal(0.1, kBlack, 0.8f));
  scene.Add(plane);

  return scene;
}

// Put boxes of varying refractive constant in the frame.
// Then move a ball across behind it - you should see the resulting time delay.
Scene dielectric_delays(int width, int height, arglist args = {}) {
  Scene scene;

  float ball_speed = 0.25f;
  if (args.size() > 0) {
    ball_speed = std::stof(args[0]);
  }

  Vec3 up(0, 1, 0);
  Vec3 right(1, 0, 0);
  Vec3 forward = Cross(up, right);
  Vec3 down = -up;
  Vec3 left = -right;
  Vec3 backward = -forward;

  float horizontal_fov = 60 * kPi / 180;
  float half_fov = horizontal_fov / 2;
  float quarter_fov = half_fov / 2;

  Camera cam(kZero3, forward, right * tanf(half_fov), width, height);
  scene.AddCamera(cam);

  scene.SetSkylight(kWhite * 0.05, (down + forward + left).NormalizedNonzero());
  scene.SetAmbientBackground(kWhite * 0.1);

  float box_height = 1.0f;
  float box_width = 1.0f;
  float box_depth = 1.0f;
  float spacing = 0.5f;

  float mid_distance = 3.0f + box_depth / 2;
  Vec3 mid_left = mid_distance * (forward + tanf(quarter_fov) * left);

  Vec3 g_center = mid_left + up * (box_height + spacing) / 2;
  Vec3 d_center = mid_left + down * (box_height + spacing) / 2;

  Material glass = kGlass;
  Material diamond = kDiamond;
  // glass.SetProperty("albedo", 0.95f);
  // diamond.SetProperty("albedo", 0.95f);

  /* Align central box axis with line to camera
  Vec3 c2r =
      (cosf(quarter_fov) * right + sinf(quarter_fov) * forward) * box_width / 2;

  Vec3 g_c2b = g_center.NormalizedNonzero() * box_depth / 2;
  Vec3 d_c2b = d_center.NormalizedNonzero() * box_depth / 2;

  Vec3 g_c2u = Cross(c2r, g_c2b).NormalizedNonzero() * box_height / 2;
  Vec3 d_c2u = Cross(c2r, d_c2b).NormalizedNonzero() * box_height / 2;

  Box glass_box(g_center, c2r, g_c2b, g_c2u, glass);
  Box diamond_box(d_center, c2r, d_c2b, d_c2u, diamond);
  */

  Vec3 c2r = right * box_width / 2;
  Vec3 c2b = forward * box_depth / 2;
  Vec3 c2u = up * box_height / 2;

  Box glass_box(g_center, c2r, c2b, c2u, glass);
  Box diamond_box(d_center, c2r, c2b, c2u, diamond);
  scene.Add(glass_box);
  scene.Add(diamond_box);

  float far_distance = mid_distance + box_depth;
  Vec3 far_left = far_distance * (forward + tanf(quarter_fov) * left);
  Vec3 far_right = far_distance * (forward + tanf(quarter_fov) * right);
  Vec3 ball_vel = ball_speed * down;
  float ball_rad = 0.2 * box_height;
  Material blue_emitter = Light(Spectrum(GaussianData(1.0f, 400, 50)));
  Sphere left_ball(far_left, ball_rad, ball_vel, 0.0f, blue_emitter);
  Sphere right_ball(far_right, ball_rad, ball_vel, 0.0f, blue_emitter);
  scene.Add(left_ball);
  scene.Add(right_ball);

  return scene;
}

// Put two emitting panels on the left & right, one red, one green.
// Have three diffuse balls move between them at different speeds.
Scene headlight_absorption(int width, int height, arglist args = {}) {
  Scene scene;

  float speed = 0.4f;
  if (args.size() > 0) {
    speed = std::stof(args[0]);
  }

  Vec3 up(0, 1, 0);
  Vec3 right(1, 0, 0);
  Vec3 forward = Cross(up, right);
  Vec3 down = -up;
  Vec3 left = -right;
  Vec3 backward = -forward;

  Vec3 cam_pos(kZero3);
  float horz_fov = 100 * kPi / 180;
  float half_fov = horz_fov / 2;
  Camera cam(cam_pos, forward, right * tanf(half_fov), width, height);
  scene.AddCamera(cam);

  // Add surrounding box, set side panels to lights
  float box_side = 5.0f;
  float half_side = box_side / 2;
  Vec3 box_O = cam_pos + forward * box_side * 0.4f;
  Box box(box_O, right * half_side,
          forward * half_side, up * half_side, Metal(0.5f, kBlack, 0.5f));

  box.right.SetMaterial(Light(Spectrum(GaussianData(1.0f, 625, 25))));
  box.left.SetMaterial(Light(Spectrum(GaussianData(1.0f, 525, 25))));
  scene.Add(box);

  // Add moving balls
  float ball_rad = box_side / 10.f;
  Material lambertian(0.95f);
  Vec3 mid_ball_pos = box_O + forward * box_side / 4;
  Vec3 top_ball_pos = mid_ball_pos + up * 3 * ball_rad;
  Vec3 bot_ball_pos = mid_ball_pos + down * 3 * ball_rad;
  Sphere mid_ball(mid_ball_pos, ball_rad, lambertian);
  float time_delay = (top_ball_pos - cam_pos).Norm() + 2*ball_rad;
  Sphere top_ball(top_ball_pos, ball_rad, speed * right, -time_delay,
                  lambertian);
  Sphere bot_ball(bot_ball_pos, ball_rad, speed * left, -time_delay,
                  lambertian);
  scene.Add(top_ball);
  scene.Add(mid_ball);
  scene.Add(bot_ball);

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
    {"headlight_absorption", &headlight_absorption}};

int main(int argc, char *argv[]) {
  // Start with standard arguments, then parse command line args.
  cmd_line_args args;
  if (argc > 1) {
    bool parsed = parse_arguments(argc, argv, args);
    if (!parsed) {
      std::cout << "Error parsing cmd line arguments. Terminating."
                << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Successfully parsed cmd line arguments." << std::endl;
  }

  // DEAL with tests if specified
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
