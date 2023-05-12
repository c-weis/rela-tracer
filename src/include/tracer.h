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

enum RenderMode { kConsole, kOutputFile, kLive };

/*
  TRACER CLASS
*/

class Tracer {
 public:
  explicit Tracer(const Scene &scene) : scene_(scene) {}

  bool RenderImage(int rays_per_pixel = 5, int depth = 3, int camera_index = 0,
                   int iterations_per_update = 1, std::string filename = "",
                   float camera_time = 0.0f) const;

  bool RenderFilm(int rays_per_pixel = 5, int depth = 3, int camera_index = 0,
                  int iterations_per_update = 1, std::string filename = "",
                  float start_time = 0.0f, float end_time = 1.0f,
                  float d_time = 0.1f) const;

 private:
  const Scene &scene_;

  rgbData TraceRay(const Line &ray, int depth, Vec3 camera_vel) const;

  float EstimateRescaleFactor(LineList *image_rays, int depth,
                              Vec3 camera_vel) const;

  void UpdateFrame(const LineList &image_rays, int depth, Vec3 cam_vel,
                   int iteration, std::vector<rgbData> &pixel_rgb) const;

  void OutputImage(int width, int height, const std::vector<rgbData> &rgb_array,
                   float rescale_factor, std::string filename) const;
};

#endif