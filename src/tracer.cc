// Copyright 2023 Christoph Weis
#include "include/tracer.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "include/materials.h"
#include "include/math.h"
#include "include/scene.h"

const int kRescaleFactorRays = 100;

bool Tracer::RenderImage(int rays_per_pixel, int depth,
                         std::vector<int> scatter_ray_counts, int camera_index,
                         RenderMode render_mode, std::string filename) const {
  // --------------------------------------------
  //    SETUP
  // --------------------------------------------
  Camera camera = scene_.GetCamera(camera_index);
  uint width = camera.GetScreenWidth();
  uint height = camera.GetScreenHeight();

  // Fill `scatter_ray_counts` up to length `depth` with 1s.
  scatter_ray_counts.insert(
      scatter_ray_counts.end(),
      std::max(depth - static_cast<int>(scatter_ray_counts.size()), 0), 1);
  std::vector<int>::const_iterator scatter_ray_count_iterator =
      scatter_ray_counts.cbegin();

  std::ofstream ofs;
  SetupOutput(width, height, render_mode, ofs, filename);

  LineList image_rays = camera.ImageRays(rays_per_pixel, r_gen_);

  float rescale_factor =
      EstimateRescaleFactor(&image_rays, kRescaleFactorRays,
                            scatter_ray_count_iterator, camera.GetVelocity());

  // Create array which we fill with the rgb values. We need this in order
  // to normalise them before converting them into RGB.
  std::vector<rgbData> rgb_array;
  rgb_array.reserve(width * height);

  // --------------------------------------------
  //    LIVE RENDER
  // --------------------------------------------

  for (auto first_pixel_ray = image_rays.cbegin();
       first_pixel_ray < image_rays.cend(); first_pixel_ray += rays_per_pixel) {
    // Average color data over rays for this pixel, turn to rgb and RGB values
    // rgb is used for rescale later, RGB is used for live output
    ColorData ray_color = kBlack;
    for (auto ray = first_pixel_ray; ray < first_pixel_ray + rays_per_pixel;
         ray++) {
      ray_color += TraceRay(*ray, depth, scatter_ray_count_iterator,
                            camera.GetVelocity());
    }
    ray_color /= static_cast<float>(rays_per_pixel);
    rgbData pixel_rgb = ray_color.To_rgb();
    rgb_array.push_back(pixel_rgb);

    RGBData pixel_RGB = pixel_rgb.ToRGB();

    OutputPixel(pixel_RGB, std::distance(image_rays.cbegin(), first_pixel_ray),
                width, height, render_mode, ofs);
  }
  ofs.close();

  // --------------------------------------------
  //    COLOR-CORRECTED RENDER
  // --------------------------------------------
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

bool Tracer::SetupOutput(int width, int height, RenderMode render_mode,
                         std::ofstream &ofs, std::string filename) const {
  switch (render_mode) {
    case kConsole:
      std::cout << "Starting console render:" << std::endl;
      break;
    case kOutputFile:
      if (filename.empty()) {
        std::cerr
            << "Error: File output mode picked, but no output file specified."
            << std::endl;
        return false;
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
      std::cerr << "Error: Specified render mode not implemented." << std::endl;
      return false;
  }
  return true;
}

float Tracer::EstimateRescaleFactor(
    LineList *image_rays, int depth,
    std::vector<int>::const_iterator scatter_ray_count_iterator,
    Vec3 camera_vel) const {
  // Estimate color rescale factor from kRescaleFactorRays random rays.
  // We assume that the average brightness is 0.5.
  std::discrete_distribution<int> random_ray_index(
      0, static_cast<int>(image_rays->size()));
  ColorData sample_color = kBlack;
  for (int sample = 0; sample < kRescaleFactorRays; sample++) {
    Line sample_ray = image_rays->at(random_ray_index(*r_gen_));
    sample_color +=
        TraceRay(sample_ray, depth, scatter_ray_count_iterator, camera_vel);
  }
  sample_color /= kRescaleFactorRays;
  rgbData sample_rgb = sample_color.To_rgb();
  float max_value =
      std::max(sample_rgb.r, std::max(sample_rgb.g, sample_rgb.b));
  if (max_value > kDivisionEpsilon) {
    return 0.5f / max_value;
  }
  return 1.0f;
}

void Tracer::OutputPixel(RGBData pixel_RGB, int pixel_index, int width,
                         int height, RenderMode render_mode,
                         std::ofstream &ofs) const {
  switch (render_mode) {
    case RenderMode::kConsole: {
      int shade = (pixel_RGB.R + pixel_RGB.G + pixel_RGB.B) / 3;
      std::cout << std::setw(3) << shade << " ";
      if ((pixel_index + 1) % width == 0) {
        std::cout << std::endl << std::flush;
      }
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

ColorData Tracer::TraceRay(const Line &ray, int depth,
                           std::vector<int>::const_iterator scatter_ray_count,
                           Vec3 camera_vel) const {
  ColorData color = ColorInStandardFrame(ray, depth, scatter_ray_count);
  color.TransformToFrame(camera_vel);
  return color;
}

ColorData Tracer::ColorInStandardFrame(
    const Line &ray, int depth,
    std::vector<int>::const_iterator scatter_ray_count) const {
  if (depth <= 0) {
    return kBlack;
  }

  OptionalHitRecord hit_ = scene_.MostRecentHit(ray);
  if (!hit_.has_value()) {
    return scene_.BackgroundColor(ray);
  }
  HitRecord hit = hit_.value();

  // extract hit data
  const Object *hit_obj = hit.obj;
  MaterialPtr mat = hit_obj->GetMaterial();
  const Line object_worldline = hit_obj->GetWorldline();
  const Vec3 object_velocity = hit_obj->GetVelocity();
  Vec4 ray_origin = hit.rf.pos.TransformedFromFrame(object_worldline);

  // de-scatter
  int scatter_rays = *scatter_ray_count;
  std::vector<ScatterData> scatter_data_obj_frame =
      mat->InverseScatter(hit.rf, scatter_rays, r_gen_);

  // recursively sum color
  ColorData color = kBlack;
  for (ScatterData s_data_of : scatter_data_obj_frame) {
    Vec3 standard_frame_vel =
        s_data_of.vel.TransformedFromFrame(object_velocity);
    Line standard_frame_ray(ray_origin, standard_frame_vel);
    ColorData pre_scatter_color =  // color of ray in object reference frame
        ColorInStandardFrame(standard_frame_ray, depth - 1,
                             scatter_ray_count + 1)
            .TransformToFrame(object_velocity);
    color += mat->ScatteredColor(pre_scatter_color, hit.rf, s_data_of.vel,
                                 s_data_of.material_data) *
             s_data_of.weight;
  }
  color /= static_cast<float>(scatter_rays);

  // add emitted color by hit object
  color += mat->EmittedColor(hit.rf);

  // transform into standard frame
  color.TransformFromFrame(hit.obj->GetVelocity());

  return color;
}
