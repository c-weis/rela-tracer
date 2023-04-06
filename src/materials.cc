// Copyright 2023 Christoph Weis
#include "include/materials.h"

#include <algorithm>
#include <functional>

#include "include/colors.h"
#include "include/math.h"

// Reflects vector in plane given the normal vector
// TODO(c): normalise `normal`?
Vec3 Reflect(const Vec3 &ray_vel, const Vec3 &normal) {
  return ray_vel - 2 * ray_vel.Dot3(normal) * normal;
}

Vec3 RandomVectorInUnitBall(std::mt19937 *r_gen) {
  Vec3 vec(1, 1, 1);
  std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
  while (vec.NormSq() > 1 || vec.NormSq() == 0.0f) {
    vec = Vec3(dist(*r_gen), dist(*r_gen), dist(*r_gen));
  }
  return vec;
}

Vec3 RandomUnitVector(std::mt19937 *r_gen) {
  Vec3 vec = RandomVectorInUnitBall(r_gen);
  return vec.NormalizedNonzero();
}

Vec3 LambertianInverseRay(Vec3 object_normal, std::mt19937 *r_gen) {
  return -(object_normal + RandomUnitVector(r_gen)).NormalizedNonzero();
}

std::vector<ScatterData> Material::InverseScatter(
    const ReferenceFrameHit &rf_hit, int scatter_rays,
    std::mt19937 *random_number_generator) const {
  std::vector<ScatterData> inverse_scatter_velocities;
  for (std::size_t component_index = 0; component_index < components_.size();
       component_index++) {
    // Scatters from each component of the material. Component weights are
    // applied here.
    WeightedComponent w_comp = components_[component_index];
    std::vector<ScatterData> inv_scatter = w_comp.comp->InverseScatter(
        rf_hit, scatter_rays, random_number_generator);
    std::for_each(inv_scatter.begin(), inv_scatter.end(),
                  [component_index, w_comp](ScatterData s_data) {
                    s_data.weight *= w_comp.weight;
                    s_data.material_data = component_index;
                  });
    inverse_scatter_velocities.insert(inverse_scatter_velocities.end(),
                                      inv_scatter.begin(), inv_scatter.end());
  }
  return inverse_scatter_velocities;
}

ColorData Material::ScatteredColor(const ColorData &pre_scatter_color,
                                   const ReferenceFrameHit &rf_hit,
                                   const Vec3 &pre_scattered,
                                   int material_data) const {
  // Material uses material_data to index the
  // scattering material responsible for each ray
  return components_[material_data].comp->ScatteredColor(pre_scatter_color,
                                                         rf_hit, pre_scattered);
}

ColorData Material::EmittedColor(const ReferenceFrameHit &rf_hit) const {
  ColorData emitted = kBlack;
  for (auto w_comp : components_) {
    emitted += w_comp.comp->EmittedColor(rf_hit) * w_comp.weight;
  }
  return emitted;
}

std::vector<ScatterData> Lambertian::InverseScatter(
    const ReferenceFrameHit &rf_hit, int scatter_rays,
    std::mt19937 *r_gen) const {
  std::vector<ScatterData> de_scatter_velocities;
  std::generate_n(
      std::back_inserter(de_scatter_velocities), scatter_rays,
      [rf_hit, r_gen]() {
        return ScatterData(LambertianInverseRay(rf_hit.normal, r_gen));
      });
  return de_scatter_velocities;
}

ColorData Lambertian::ScatteredColor(const ColorData &pre_scatter_color,
                                     const ReferenceFrameHit &rf_hit,
                                     const Vec3 &pre_scattered,
                                     int material_data) const {
  return pre_scatter_color * albedo_;
}

std::vector<ScatterData> Dielectric::InverseScatter(
    const ReferenceFrameHit &rf_hit, int scatter_rays,
    std::mt19937 *random_number_generator) const {
  // Adapted from “Ray Tracing in One Weekend.”
  // raytracing.github.io/books/RayTracingInOneWeekend.html
  // (we need to keep track of the speed of the ray.)
  Vec3 ray_vel = rf_hit.scattered;
  Vec3 normal = rf_hit.normal;

  float n_ratio;
  if (rf_hit.outside) {
    n_ratio = 1 / n_;
  } else {
    n_ratio = n_;
  }
  float n_ratio_2 = n_ratio * n_ratio;

  float cos_theta = ray_vel.NormalizedNonzero().Dot3(normal);
  float new_sin_theta = sqrtf(1 - cos_theta * cos_theta) * n_ratio;

  if (new_sin_theta > 1.0f) {
    // reflect fully
    return {ScatterData(Reflect(ray_vel, normal),
                        static_cast<float>(scatter_rays))};
  }

  // Schlick's approximation
  float r0 = (1 - n_ratio) / (1 + n_ratio);
  r0 = r0 * r0;
  float reflectance = r0 + (1 - r0) * powf(1 - cos_theta, 5);

  ScatterData reflected(Reflect(ray_vel, normal), scatter_rays * reflectance);

  Vec3 parallel = Dot3(ray_vel, normal) * normal;
  Vec3 new_perpendicular = n_ratio_2 * (ray_vel - parallel);
  Vec3 new_parallel =
      sqrtf(ray_vel.NormSq() * n_ratio_2 - new_perpendicular.NormSq()) * normal;
  Vec3 refracted_vel = new_perpendicular + new_parallel;

  ScatterData refracted(refracted_vel, scatter_rays * (1 - reflectance));
  return {reflected, refracted};
}

ColorData Dielectric::ScatteredColor(const ColorData &pre_scatter_color,
                                     const ReferenceFrameHit &rf_hit,
                                     const Vec3 &pre_scattered,
                                     int material_data) const {
  return pre_scatter_color;
}

ColorData MonochromaticLight::EmittedColor(
    const ReferenceFrameHit &rf_hit) const {
  Vec4 k(omega_, omega_ * rf_hit.scattered / rf_hit.scattered.NormSq());
  return ColorData(
      {LightRay(brightness_ * Dot3(rf_hit.scattered, rf_hit.normal), k)});
}

std::vector<ScatterData> Metal::InverseScatter(const ReferenceFrameHit &rf_hit,
                                               int scatter_rays,
                                               std::mt19937 *r_gen) const {
  // Following “Ray Tracing in One Weekend.”
  // raytracing.github.io/books/RayTracingInOneWeekend.html
  Vec3 reflected_ray = Reflect(rf_hit.scattered, rf_hit.normal);
  if (fuzz_ == 0.0f) {
    return {ScatterData(reflected_ray, static_cast<float>(scatter_rays))};
  }
  std::vector<ScatterData> prescattered_rays;
  prescattered_rays.assign(scatter_rays, ScatterData(reflected_ray));
  for (ScatterData &w_vel : prescattered_rays) {
    w_vel.vel = w_vel.vel + fuzz_ * RandomVectorInUnitBall(r_gen);
  }
  return prescattered_rays;
}

ColorData Metal::ScatteredColor(const ColorData &pre_scatter_color,
                                const ReferenceFrameHit &rf_hit,
                                const Vec3 &pre_scattered,
                                int material_data) const {
  return pre_scatter_color * albedo_;
}
