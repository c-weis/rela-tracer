// Copyright 2023 Christoph Weis
#include "include/tracer.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "include/materials.h"
#include "include/math.h"
#include "include/scene.h"

bool Tracer::RenderImage(int rays_per_pixel, int depth,
                         std::vector<int> scatter_ray_counts, int camera_index,
                         RenderMode render_mode, std::string filename) const {
  uint width = scene_.GetCamera(camera_index).ScreenWidth();
  uint height = scene_.GetCamera(camera_index).ScreenHeight();

  // Fill `scatter_ray_counts` up to length `depth` with 1s.
  scatter_ray_counts.insert(
      scatter_ray_counts.end(),
      std::max(depth - static_cast<int>(scatter_ray_counts.size()), 0), 1);
  std::vector<int>::const_iterator scatter_ray_count_iterator =
      scatter_ray_counts.cbegin();

  // Create array which we fill with the rgb values. We need this in order
  // to normalise them before converting them into RGB.
  std::vector<rgbData> rgb_array;
  rgb_array.reserve(width * height);

  std::ofstream ofs;
  switch (render_mode) {
    case kConsole:
      std::cout << "Starting console render:" << std::endl;
      break;
    case kOutputFile:
      if (filename.empty()) {
        std::cerr
            << "Error: File output mode picked, but no output file specified."
            << std::endl;
        return EXIT_FAILURE;
      }
      std::cout << "Outputting to file: " << filename << "_live.ppm"
                << std::endl;
      ofs.open(filename + "_live.ppm",
               std::ios_base::out | std::ios_base::binary);
      ofs << "P6" << std::endl
          << width << ' ' << height << std::endl
          << "255" << std::endl;
      break;
    default:
      std::cerr << "Error: Specified render mode not implemented.";
      return EXIT_FAILURE;
  }

  LineList image_rays =
      scene_.GetCamera(camera_index).ImageRays(rays_per_pixel, r_gen_);

  // Iterate over pixels, perform live render.
  for (auto first_pixel_ray = image_rays.cbegin();
       first_pixel_ray < image_rays.cend(); first_pixel_ray += rays_per_pixel) {
    // Average color data over rays for this pixel, turn to rgb and RGB values
    // rgb is used for rescale later, RGB is used for live output
    ColorData ray_color = kBlack;
    for (auto ray = first_pixel_ray; ray < first_pixel_ray + rays_per_pixel;
         ray++) {
      ray_color += TraceRay(*ray, depth, scatter_ray_count_iterator);
    }
    ray_color /= static_cast<float>(rays_per_pixel);
    rgbData pixel_rgb = ray_color.To_rgb();
    rgb_array.push_back(pixel_rgb);
    // TODO(c): estimate good rescale factor (e.g. by sampling random pixels
    //          at the beginning). Use it to get better live render.
    RGBData pixel_RGB = pixel_rgb.ToRGB();

    switch (render_mode) {
      case RenderMode::kConsole: {
        int shade = (pixel_RGB.R + pixel_RGB.G + pixel_RGB.B) / 3;
        std::cout << std::setw(3) << shade << " ";
        /*
        if (shade >= 128) {
          std::cout << 0;
        } else {
          std::cout << " ";
        }
        */
        if ((std::distance(image_rays.cbegin(), first_pixel_ray) /
                 rays_per_pixel +
             1) %
                width ==
            0)
          std::cout << std::endl << std::flush;
        break;
      }
      case RenderMode::kOutputFile: {
        ofs << static_cast<char>(pixel_RGB.R) << static_cast<char>(pixel_RGB.G)
            << static_cast<char>(pixel_RGB.B);
        break;
      }
      default:
        std::cout << ".";
        break;
    }
  }
  ofs.close();

  // TEST
  std::vector<float> min_rgb_array;
  std::transform(rgb_array.cbegin(), rgb_array.cend(),
                 std::back_inserter(min_rgb_array), [](const rgbData rgb) {
                   return std::min(std::min(rgb.r, rgb.g), rgb.b);
                 });
  float min_rgb =
      *(std::min_element(min_rgb_array.cbegin(), min_rgb_array.cend()));
  if (render_mode == RenderMode::kOutputFile) {
    // Perform normalised render:
    // We find the maximum rgb value, compute normalised RGB and output the
    // result.
    std::vector<float> max_rgb_array;
    std::transform(rgb_array.cbegin(), rgb_array.cend(),
                   std::back_inserter(max_rgb_array), [](const rgbData rgb) {
                     return std::max(std::max(rgb.r, rgb.g), rgb.b);
                   });
    float max_rgb =
        *(std::max_element(max_rgb_array.cbegin(), max_rgb_array.cend()));
    const float rescale_factor = 255.0f / max_rgb;
    std::vector<RGBData> RGB_array;
    RGB_array.reserve(width * height);
    std::transform(rgb_array.cbegin(), rgb_array.cend(),
                   std::back_inserter(RGB_array),
                   [rescale_factor](const rgbData rgb) {
                     return rgb.ToRGB(rescale_factor);
                   });

    std::cout << "Outputting to file: " << filename << ".ppm" << std::endl;
    ofs.open(filename + ".ppm", std::ios_base::out | std::ios_base::binary);
    ofs << "P6" << std::endl
        << width << ' ' << height << std::endl
        << "255" << std::endl;

    for (auto [R, G, B] : RGB_array) {
      ofs << static_cast<char>(R) << static_cast<char>(G)
          << static_cast<char>(B);
    }
    ofs.close();
  }

  return EXIT_SUCCESS;
}

std::vector<Vec3> Tracer::InverseScatter(HitRecord hit,
                                         int scattered_ray_count) const {
  return hit.obj->mat_->InverseScatter(hit.rf.scattered, hit.rf.normal,
                                       scattered_ray_count, r_gen_);
}

ColorData Tracer::ColorInStandardFrame(ScatterNode root) const {
  // Average children's colours recursively.
  // Apply color transforms.

  const Line worldline = root.hit.obj->worldline_;
  const Material *mat = root.hit.obj->mat_;
  const ReferenceFrameHit rf_hit = root.hit.rf;

  // average children
  ColorData color = kBlack;
  for (ScatterNode child_node : root.children) {
    // obtain color from child (recursively) &
    // apply Lorentz boost to rest frame of hit object
    ColorData pre_scattered_color =
        ColorInStandardFrame(child_node).TransformToFrame(worldline.vel);
    color += mat->ScatteredColor(rf_hit, child_node.ray_vel_in_parent_rf,
                                 pre_scattered_color) /
             root.children.size();
  }

  color += mat->EmittedColor(rf_hit);

  // de-boost to standard frame
  return color.TransformFromFrame(worldline.vel);
}

OptionalScatterNode Tracer::RecursiveTraceRay(
    const Line &ray, int depth,
    std::vector<int>::const_iterator scatter_ray_count,
    Line parent_restframe_worldline) const {
  if (depth <= 0) {
    return std::nullopt;
  }

  // find nearest intersection
  OptionalHitRecord hit_ = scene_.MostRecentHit(ray);

  if (!hit_) {
    return std::nullopt;
  }

  HitRecord hit = hit_.value();

  /*
    TODO(c): make scatter_rays an argument in function, allow scheduling
  */
  int scatter_rays = *scatter_ray_count;

  // de-scatter
  ScatterNode scatter_node{.ray_vel = ray.vel,
                           .ray_vel_in_parent_rf = ray.vel.TransformedToFrame(
                               parent_restframe_worldline.vel),
                           .hit = hit};

  std::vector<Vec3> obj_frame_pre_scatter_velocities =
      InverseScatter(hit, scatter_rays);

  Line restframe_worldline = hit.obj->worldline_;
  Vec4 ray_origin = hit.rf.pos.TransformedFromFrame(hit.obj->worldline_);

  for (Vec3 obj_frame_vel : obj_frame_pre_scatter_velocities) {
    Vec3 standard_frame_vel =
        obj_frame_vel.TransformedFromFrame(restframe_worldline.vel);
    Line ray(ray_origin, standard_frame_vel);
    OptionalScatterNode child_node = RecursiveTraceRay(
        ray, depth - 1, scatter_ray_count + 1, restframe_worldline);
    if (child_node.has_value()) {
      scatter_node.children.push_back(child_node.value());
    }
  }

  return scatter_node;
}

ColorData Tracer::TraceRay(
    const Line &ray, int depth,
    std::vector<int>::const_iterator scatter_ray_count) const {
  Line standard_worldine = Line(kZero4, kZero3);
  OptionalScatterNode scatter_tree_ =
      RecursiveTraceRay(ray, depth, scatter_ray_count, standard_worldine);

  if (!scatter_tree_.has_value()) {
    // TEST: depth of scatter tree
    // std::cout << "-1 ";
    return kBlack;
  }
  ScatterNode scatter_tree = scatter_tree_.value();

  // TEST: get depth of scatter tree
  // std::cout << scatter_depth(scatter_tree) << " ";
  // if (scatter_depth(scatter_tree) > 0) {
  //   std::cout << "(!)"
  //             << scatter_tree.hit.rf.pos -
  //             scatter_tree.children[0].hit.rf.pos
  //             << " "
  //             << Dot3(scatter_tree.children[0].hit.rf.normal,
  //                     scatter_tree.children[0].hit.rf.scattered)
  //             << std::endl;
  // }

  ColorData color = ColorInStandardFrame(scatter_tree);

  // TODO(c): implement moving cam here by boosting to cam rf

  return color;
}

int scatter_depth(const ScatterNode &root, int start_depth) {
  int max_depth = start_depth;
  for (ScatterNode child : root.children) {
    max_depth = std::max(max_depth, scatter_depth(child, start_depth + 1));
  }
  return max_depth;
}
