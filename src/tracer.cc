// Copyright 2023 Christoph Weis
#include "include/tracer.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "include/bitmap_image.h"
#include "include/materials.h"
#include "include/math.h"
#include "include/scene.h"
#include "tracer.h"

const int kRescaleFactorRays = 1000;

// Updates the given pixel color array by tracing along the image rays with
// the specified parameters.
// Note that iteration=0 automatically overwrites the rgb_array, which allows
// easy reuse of the same array for independent frames in a film.
void Tracer::UpdateFrame(const LineList &image_rays, int depth, Vec3 cam_vel,
                         int iteration, std::vector<rgbData> &rgb_array) const {
  for (int pixel = 0; pixel < image_rays.size(); pixel++) {
    Line ray = image_rays[pixel];

    rgbData pixel_rgb = TraceRay(ray, depth, cam_vel);

    rgb_array[pixel] = (rgb_array[pixel] * iteration + pixel_rgb) /
                       (static_cast<float>(iteration) + 1.0f);
  }
}

std::string AssembleFilename(std::string prefix, std::string name, int number,
                             int pad_to_digits, std::string postfix) {
  std::string framenumber = std::to_string(number);
  framenumber =
      std::string(pad_to_digits - framenumber.size(), '0') + framenumber;
  return prefix + name + "_" + framenumber + postfix;
}

float MaxValue(const std::vector<rgbData> &rgb_array) {
  std::vector<float> max_rgb_array;
  std::transform(rgb_array.cbegin(), rgb_array.cend(),
                 std::back_inserter(max_rgb_array), [](const rgbData rgb) {
                   return std::max(std::max(rgb.r, rgb.g), rgb.b);
                 });
  return *(std::max_element(max_rgb_array.cbegin(), max_rgb_array.cend()));
}

// Verbosely renders an image with the given parameters.
bool Tracer::RenderImage(int rays_per_pixel, int depth, int camera_index,
                         int iterations_per_update, std::string filename,
                         float camera_time, float rescale_factor,
                         bool output_normalised_render) const {
  // --------------------------------------------
  //    SETUP
  // --------------------------------------------
  Camera camera = scene_.GetCamera(camera_index);
  unsigned int width = camera.GetScreenWidth();
  unsigned int height = camera.GetScreenHeight();

  LineList image_rays = camera.ImageRays(camera_time);

  bool estimate_rescale_factor = rescale_factor < 0.0f;
  if (estimate_rescale_factor) {
    rescale_factor = EstimateRescaleFactor(
        &image_rays, depth, camera.GetVelocity(), kRescaleFactorRays);
  }

  // Create array which we fill with the rgb values. We need this in order
  // to normalise them before converting them into RGB.
  std::vector<rgbData> rgb_array;
  rgb_array.assign(width * height, rgbData(0, 0, 0));

  // Complete file names
  std::string live_filename = "images/" + filename + ".bmp";
  std::string normalised_filename = "images/n_" + filename + ".bmp";

  // --------------------------------------------
  //    LIVE RENDER
  // --------------------------------------------
  std::cout << "Rendering scene with:" << std::endl;
  std::cout << "depth: " << depth << ", rays per pixel:" << rays_per_pixel
            << std::endl;
  std::cout << "Rescale factor: " << rescale_factor;
  if (estimate_rescale_factor) {
    std::cout << " (estimated)";
  }
  std::cout << std::endl;

  std::cout << "Output to file " << live_filename << " every "
            << iterations_per_update << " iterations." << std::endl;
  std::cout << "Starting render." << std::endl;
  for (int iteration = 0; iteration < rays_per_pixel; iteration++) {
    // refresh image rays
    image_rays = camera.ImageRays(camera_time);

    UpdateFrame(image_rays, depth, camera.GetVelocity(), iteration, rgb_array);

    if ((iteration + 1) % iterations_per_update == 0) {
      OutputImage(width, height, rgb_array, rescale_factor, live_filename);
      std::cout << "Done " << (iteration + 1) << "/" << rays_per_pixel
                << std::endl;
    }
  }

  if (!output_normalised_render) {
    return EXIT_SUCCESS;
  }
  // --------------------------------------------
  //    COLOR-CORRECTED RENDER
  // --------------------------------------------
  // Perform normalised render:
  // We find the maximum rgb value, compute normalised RGB and output.
  std::cout << "Normalising render, output to " << normalised_filename << "."
            << std::endl;
  rescale_factor = 255.0f / MaxValue(rgb_array);

  OutputImage(width, height, rgb_array, rescale_factor, normalised_filename);
  return EXIT_SUCCESS;
}

/*
  TODO(c): add more documentation
*/
bool Tracer::RenderFilm(int rays_per_pixel, int depth, int camera_index,
                        int iterations_per_update, std::string filename,
                        float start_time, float end_time, float d_time,
                        int start_frame, int end_frame, bool preview_only,
                        float rescale_factor,
                        bool output_normalised_renders) const {
  // --------------------------------------------
  //    SETUP
  // --------------------------------------------
  Camera camera = scene_.GetCamera(camera_index);
  unsigned int width = camera.GetScreenWidth();
  unsigned int height = camera.GetScreenHeight();

  const int frames = (end_time - start_time) / d_time + 1;

  int last_frame = frames;
  if (end_frame != -1) {
    last_frame = end_frame;
  }

  bool estimate_rescale_factor = rescale_factor < 0.0f;
  if (estimate_rescale_factor) {
    // Estimate rescale factor
    rescale_factor = 0.0f;
    for (int estimator_frame = start_frame; estimator_frame < last_frame;
         estimator_frame++) {
      float estimator_time = start_time + d_time * estimator_frame;
      LineList image_rays = camera.ImageRays(estimator_time);
      rescale_factor += EstimateRescaleFactor(
          &image_rays, depth, camera.GetVelocity(), kRescaleFactorRays);
    }
    rescale_factor /= last_frame - start_frame;
    std::cout << "Estimated rescale factor: " << rescale_factor << std::endl;
  }

  // --------------------------------------------
  //    PREVIEW STAGE
  // --------------------------------------------
  // Quickly render a preview of each frame
  const int preview_rays_per_pixel = 1;
  const int preview_depth = 5;

  const int frame_digits = static_cast<int>(std::floor(std::log10(frames))) + 1;

  std::cout << "Outputting preview images to vidframes/" << filename << "_XXX."
            << std::endl;

  // Create array which we fill with the rgb values. We need this in order
  // to normalise them before converting them into RGB.
  std::vector<rgbData> rgb_array;
  rgb_array.assign(width * height, rgbData(0, 0, 0));

  // Loop through preview frames
  for (int preview_frame = start_frame; preview_frame < last_frame;
       preview_frame++) {
    float preview_time = start_time + d_time * preview_frame;
    std::string preview_filename = AssembleFilename(
        "vidframes/", filename, preview_frame, frame_digits, ".bmp");
    for (int iteration = 0; iteration < preview_rays_per_pixel; iteration++) {
      LineList image_rays = camera.ImageRays(preview_time);
      UpdateFrame(image_rays, preview_depth, camera.GetVelocity(), iteration,
                  rgb_array);
    }
    OutputImage(width, height, rgb_array, rescale_factor, preview_filename);
    std::cout << "Finished preview frame " << preview_frame + 1 << "."
              << std::endl;
  }
  std::cout << "---------------------------------" << std::endl;

  if (preview_only) {
    std::cout << "Program was run with preview_only. Done." << std::endl;
    return EXIT_SUCCESS;
  }

  // --------------------------------------------
  //    LIVE RENDER
  // --------------------------------------------
  std::cout << "Starting render." << std::endl;
  std::cout << "Rendering scene with:" << std::endl;
  std::cout << "depth: " << depth << ", rays per pixel:" << rays_per_pixel
            << std::endl;

  std::cout << "Rescale factor: " << rescale_factor;
  if (estimate_rescale_factor) {
    std::cout << " (estimated)";
  }
  std::cout << std::endl;

  std::cout << "Updating files every " << iterations_per_update
            << " iterations." << std::endl;
  for (int frame = start_frame; frame < last_frame; frame++) {
    float camera_time = start_time + d_time * frame;
    std::string complete_filename =
        AssembleFilename("vidframes/", filename, frame, frame_digits, ".bmp");

    for (int iteration = 0; iteration < rays_per_pixel; iteration++) {
      // refresh image rays
      LineList image_rays = camera.ImageRays(camera_time);

      UpdateFrame(image_rays, depth, camera.GetVelocity(), iteration,
                  rgb_array);

      if ((iteration + 1) % iterations_per_update == 0 ||
          iteration + 1 == rays_per_pixel) {
        OutputImage(width, height, rgb_array, rescale_factor,
                    complete_filename);
        std::cout << "Frame " << frame << " iteration " << (iteration + 1)
                  << "/" << rays_per_pixel << std::endl;
      }
    }

    if (output_normalised_renders) {
      // --------------------------------------------
      //    COLOR-CORRECTED RENDER
      // --------------------------------------------
      // Perform normalised render:
      // Find the maximum rgb value, compute normalised RGB and output.
      std::cout << "Normalising frame " << frame << ".";
      float frame_rescale_factor = 255.0f / MaxValue(rgb_array);
      OutputImage(width, height, rgb_array, frame_rescale_factor,
                  AssembleFilename("vidframes/", "n_" + filename, frame,
                                   frame_digits, ".bmp"));
      std::cout << " Done." << std::endl
                << "---------------------------------" << std::endl;
    }
  }

  return EXIT_SUCCESS;
}

float Tracer::EstimateRescaleFactor(LineList *image_rays, int depth,
                                    Vec3 camera_vel,
                                    int estimation_rays) const {
  // Estimate color rescale factor from kRescaleFactorRays random rays.
  // We assume that the average brightness is 0.5.
  rgbData sample_color;
  for (int sample = 0; sample < estimation_rays; sample++) {
    Line sample_ray =
        image_rays->at(static_cast<int>(RandomReal() * image_rays->size()));
    sample_color += TraceRay(sample_ray, depth, camera_vel);
  }
  sample_color /= estimation_rays;
  float max_value =
      std::max(sample_color.r, std::max(sample_color.g, sample_color.b));
  if (max_value > kDivisionEpsilon) {
    return 0.5f / max_value * 255;
  }
  return 255.0f;
}

void Tracer::OutputImage(int width, int height,
                         const std::vector<rgbData> &rgb_array,
                         float rescale_factor, std::string filename) const {
  BitmapImage img(width, height);

  for (unsigned int y = 0; y < height; y++) {
    for (unsigned int x = 0; x < width; x++) {
      RGBData pixel = rgb_array[y * width + x].ToRGB(rescale_factor);
      if (pixel.R + pixel.G + pixel.B > 0) {
        // printf("pixel %d %d %d @ %d %d\n", pixel.R, pixel.G, pixel.B, x, y);
      }
      img.set_pixel(x, y, pixel.R, pixel.G, pixel.B);
    }
  }

  img.save_image(filename);
}

rgbData Tracer::TraceRay(const Line &ray_0, int depth, Vec3 camera_vel) const {
  rgbData color;
  Line ray(ray_0);
  SpectrumTransform cumulative_transform;
  cumulative_transform.ApplyTransformationToFrame(ray.vel, camera_vel);

  for (int d = 0; d < depth; d++) {
    OptionalHitRecord hit_ = scene_.MostRecentHit(ray);
    if (!hit_.has_value()) {  // end of ray
      return color +
             cumulative_transform.ColorFrom(scene_.EscapedRayColor(ray));
    }
    HitRecord hit = hit_.value();

    // extract hit data
    const Object *hit_obj = hit.obj;
    const Material *mat = hit_obj->GetMaterial();
    const Line object_worldline = hit_obj->GetWorldline();
    const Vec3 object_velocity = hit_obj->GetVelocity();
    Vec4 ray_origin = hit.rf.pos.TransformedFromFrame(object_worldline);

    // DE-SCATTER:
    // adjust cumulative_transform, add color & absorption, compute new ray
    cumulative_transform.ApplyTransformationFromFrame(hit.rf.scattered,
                                                      object_velocity);
    // apply color data
    color += cumulative_transform.ColorFrom(mat->EmissionSpectrum(hit.rf));
    cumulative_transform.ApplyAbsorption(mat->AbsorptionCurve(hit.rf));

    // de-scatter
    ScatterData scatter_data =
        mat->InverseScatter(hit.rf, cumulative_transform);
    Vec3 ray_vel = scatter_data.vel_.VelTransformedFromFrame(object_velocity);
    ray = Line(ray_origin, ray_vel);
    cumulative_transform.ApplyTransformationToFrame(ray_vel, object_velocity);
  }

  return color;
}
