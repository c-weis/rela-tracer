// Copyright 2023 Christoph Weis
#pragma once

#include <array>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "include/materials.h"
#include "include/math.h"
#include "include/scene.h"

enum RenderMode { kConsole, kOutputFile, kLive };

/*
  SCATTER DATA
*/

// Scatter data naturally assembles into a rooted tree:
// We have the final ray hitting the camera - the root node.
// At each scatter point, the pre-scattered rays hit the
// point of scattering and combine into the scattered ray.
// (We traverse the scattter tree backwards in time.)
struct ScatterNode {
  // Outgoing ray velocity in standard frame
  Vec3 ray_vel;
  // Outgoing ray velocity in restframe of parent
  Vec3 ray_vel_in_parent_rf;
  // Hit object, position, normal, scattered ray in obj restframe
  HitRecord hit;
  // Children (ie. pre-scattered rays)
  std::vector<ScatterNode> children;
};
typedef std::optional<ScatterNode> OptionalScatterNode;

// Return maximum depth of scatter tree
int scatter_depth(const ScatterNode &root, int start_depth = 0);

/*
  TRACER CLASS
*/

class Tracer {
 public:
  explicit Tracer(const Scene &scene)
      : scene_(scene), r_gen_(new std::mt19937(42)) {}

  bool RenderImage(int rays_per_pixel = 5, int depth = 3, std::vector<int> scatter_ray_counts = {5}, int camera_index = 0,
                   RenderMode render_mode = RenderMode::kConsole,
                   std::string filename = "") const;

 private:
  const Scene &scene_;

  std::mt19937 *r_gen_;

  ColorData TraceRay(const Line &r, int depth, std::vector<int>::const_iterator scatter_ray_count) const;

  ColorData ColorInStandardFrame(ScatterNode) const;

  OptionalScatterNode RecursiveTraceRay(const Line &ray, int depth,
                                        std::vector<int>::const_iterator scatter_ray_count,
                                        Line parent_restframe_worldlline) const;

  std::vector<Vec3> InverseScatter(HitRecord hit, int scattered_rays) const;
};
