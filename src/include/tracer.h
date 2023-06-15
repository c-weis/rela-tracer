// Copyright 2023 Christoph Weis
#ifndef TRACER_H_INCLUDED
#define TRACER_H_INCLUDED

#include <array>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "colors.h"
#include "materials.h"
#include "math.h"
#include "scene.h"

/*
  TRACER CLASS

  Takes in a scene of uniformly moving objects and camera data.
  Produces a render by the usual ray-tracing procedure,
  modified to incorporate special relativity.

  Tracer can output a single image or a series of images for a film.

  Use this class as follows:
  Define a scene, adding objects and camera(s):
    Scene scene;
    ADD OBJECTS AND CAMERA
  Instantiate tracer:
    Tracer tracey(scene);
  Run tracer, outputting a single image or frames of a film:
    tracey.RenderImage();
    OR
    tracey.RenderFilm();
*/
class Tracer {
 public:
  // Constructor: needs scene complete with objects and camera(s).
  explicit Tracer(const Scene &scene) : scene_(scene) {}

  // Render a single frame and output to "images/"+`filename`+".bmp".
  // Args:
  //  rays_per_pixel: number of rays to average over for each pixel
  //  depth: number of bounces to follow for each ray
  //  camera_index: Tracer allows multiple cameras per scene, this specifies
  //    which camera to use (0-indexed)
  //  iterations_per_update: how many rays per pixel are traced between updates
  //    of the image file
  //  filename: specifies the name of the output file
  //  camera_time: time the picture is taken, given in the camera restframe
  //  rescale_factor: Overall factor multiplying RGB values before output.
  //    Setting this to -1 will trigger a rescale factor estimation.
  //  output_normalised_render: setting this true will trigger a "normalised
  //    render", an output where the rescale factor is adjusted to avoid any
  //    clipping of color outputs to a file "n_"+`filename` + ".bmp".
  bool RenderImage(int rays_per_pixel = 5, int depth = 3, int camera_index = 0,
                   int iterations_per_update = 1, std::string filename = "",
                   float camera_time = 0.0f, float rescale_factor = -1.0f,
                   bool output_normalised_render = true) const;

  // Render frames of a film, output to "vidframes/"+`filename`+frame+".bmp"
  // Args:
  //  rays_per_pixel: number of rays to average over for each pixel
  //  depth: number of bounces to follow for each ray
  //  camera_index: Tracer allows multiple cameras per scene, this specifies
  //    which camera to use (0-indexed)
  //  iterations_per_update: how many rays per pixel are traced between updates
  //    of each image file
  //  filename: specifies the name of the output files
  //  start_time: time the first frame is taken, in the camera restframe
  //  end_time: time the last frame is taken, in the camera restframe
  //  d_time: time between frames, in the camera restframe
  //  start_frame: index of first frame to be rendered (useful for parallel
  //    rendering)
  //  end_frame: index of last frame to be rendered (set to -1 to
  //    render all following frames)
  //  preview_only: if true, will only trace one ray per pixel, with depth 5
  //  rescale_factor: Overall factor multiplying RGB values before output.
  //    Setting this to -1 will trigger a rescale factor estimation.
  //  output_normalised_render: setting this true will trigger a "normalised
  //    render", an output where the rescale factor is adjusted to avoid any
  //    clipping of color outputs to files "n_"+`filename` + frame + ".bmp".
  bool RenderFilm(int rays_per_pixel = 5, int depth = 3, int camera_index = 0,
                  int iterations_per_update = 1, std::string filename = "",
                  float start_time = 0.0f, float end_time = 1.0f,
                  float d_time = 0.1f, int start_frame = 0, int end_frame = -1,
                  bool preview_only = false, float rescale_factor = -1.0f,
                  bool output_normalised_renders = false) const;

 private:
  // The scene to be rendered.
  const Scene &scene_;

  // Estimates a value for `rescale_factor`, which multiplies RGB values before
  // output to bitmap.
  float EstimateRescaleFactor(LineList *image_rays, int depth, Vec3 camera_vel,
                              int estimation_rays) const;

  // Returns color picked up by a ray.
  // Args:
  //  ray: the incoming light ray, specified in the standard frame
  //  depth: maximum bounces to trace the ray for
  //  camera_vel: velocity of camera in standard frame
  rgbData TraceRay(const Line &ray, int depth, Vec3 camera_vel) const;

  // Updates `pixel_rgb` to current iteration.
  // If iteration = 0, this amounts to replacing `pixel_rgb` with the rgb values
  // associated to `image_rays`, `depth` and `cam_vel`, otherwise it is the
  // average of all `pixel_rgb` values collected so far.
  void UpdateFrame(const LineList &image_rays, int depth, Vec3 cam_vel,
                   int iteration, std::vector<rgbData> *pixel_rgb) const;

  // Outputs `rgb_array` to bitmap file.
  void OutputImage(int width, int height, const std::vector<rgbData> &rgb_array,
                   float rescale_factor, std::string filename) const;
};

#endif
