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
  explicit Tracer(const Scene &scene)
      : scene_(scene), r_gen_(new std::mt19937(42)) {}

  bool RenderImage(int rays_per_pixel = 5, int depth = 3,
                   std::vector<int> scatter_ray_counts = {5},
                   int camera_index = 0,
                   RenderMode render_mode = RenderMode::kConsole,
                   std::string filename = "") const;

 private:
  const Scene &scene_;

  std::mt19937 *r_gen_;

  ColorData ColorInStandardFrame(
      const Line &ray, int depth,
      std::vector<int>::const_iterator scatter_ray_count) const;

  ColorData TraceRay(const Line &r, int depth,
                     std::vector<int>::const_iterator scatter_ray_count,
                     Vec3 camera_vel) const;

  std::vector<ScatterData> InverseScatter(HitRecord hit, int scattered_rays) const;

  bool SetupOutput(int width, int height, RenderMode render_mode,
                   std::ofstream &ofs, std::string filename) const;

  float EstimateRescaleFactor(
      LineList *image_rays, int depth,
      std::vector<int>::const_iterator scatter_ray_count_iterator,
      Vec3 camera_vel) const;

  void OutputPixel(RGBData pixel_RGB, int pixel_index, int width, int height,
                   RenderMode render_mode, std::ofstream &ofs) const;
};

#endif