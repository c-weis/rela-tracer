// Copyright 2023 Christoph Weis
#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// List of arguments
typedef std::vector<std::string> arglist;

// Command line argument parser.
class cmd_line_args {
 public:
  // Scene and output args
  std::string scene_name = "headlight_absorption";  // name of example scene
  std::string filename = "output";                  // name of output filename
  bool film = false;  // set true if rendering frames of a film
  bool preview_only =
      false;          // set true if only preview picture(s) shall be generated
  bool test = false;  // run tests instead of render

  // Standard rendering arguments
  int width = 500;                // image width in pixels
  int height = 500;               // image height in pixels
  int rays_per_pixel = 100;       // nr of rays to average over per pixel
  int depth = 8;                  // max number of bounces to trace rays
  int iterations_per_update = 1;  // number of passes between image updates
  int camera_index = 0;  // index of camera (a scene may contain several)

  // Time arguments
  float time = 0.0f;  // time of picture in (std frame), start_time for film
  float end_time =
      1.0f;  // time of last frame for film (std frame), ignored for image
  float d_time =
      0.1f;  // time between frames for film (std frame) ignored for image

  // Parallel processing helps
  int start_frame = 0;        // first frame to render
  int end_frame = -1;         // last frame to render (-1 = all remaining)
  float rescale_factor = -1;  // multiplier for RGB values (-1 = estimate it)

  // Scene arguments
  arglist scene_args = {};  // scene-specific arguments

  // Parses command line arguments
  bool parse_arguments(int argc, char *argv[]);

 private:
  // We split up arguments by type and give names and aliases.

  // Integer arguments
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

  // Float arguments
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

  // String arguments
  std::map<std::string, std::string &> string_args = {
      {"--scene_name", scene_name},
      {"-sn", scene_name},
      {"--filename", filename},
      {"-fn", filename}};

  // Bool arguments
  std::map<std::string, bool &> bool_args = {{"--test", test},
                                             {"-t", test},
                                             {"--film", film},
                                             {"-f", film},
                                             {"--preview_only", preview_only},
                                             {"-po", preview_only}};

  // Triggers to start list of scene-specific arguments
  std::vector<std::string> scene_args_triggers = {"--scene_args", "--args",
                                                  "-sa"};
};

// Parses command line arguments
bool cmd_line_args::parse_arguments(int argc, char *argv[]) {
  // skip first argument: it's the filename
  for (int i = 1; i < argc; i++) {
    // check if it's an integer argument
    {
      auto int_iter = int_args.find(argv[i]);
      if (int_iter != int_args.end()) {
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
      auto float_iter = float_args.find(argv[i]);
      if (float_iter != float_args.end()) {
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
      auto string_iter = string_args.find(argv[i]);
      if (string_iter != string_args.end()) {
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
      auto bool_iter = bool_args.find(argv[i]);
      if (bool_iter != bool_args.end()) {
        bool_iter->second = true;
        continue;
      }
    }

    // check if it's a trigger to start reading scene-specific arguments
    {
      auto sa_trigger_iter = std::find(scene_args_triggers.cbegin(),
                                       scene_args_triggers.cend(), argv[i]);
      if (sa_trigger_iter < scene_args_triggers.cend()) {
        // add all following arguments into the scene_args vector
        for (i = i + 1; i < argc; i++) {
          scene_args.push_back(argv[i]);
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

#endif
