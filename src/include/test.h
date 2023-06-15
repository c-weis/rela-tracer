// Copyright 2023 Christoph Weis
#ifndef TEST_H_INCLUDED
#define TEST_H_INCLUDED

#include "tracer.h"

// All-round testing class for rela-tracer components.
class Tester {
 public:
  // Test random number and vector generators
  bool TestRandomness();
  // Output color test pictures
  bool TestColors();
  // Test Lorentz transformations of positions and velocities
  bool TestLineLorentzTransforms();
  // Test Lorentz transformation of color and brightness
  bool TestColorLorentzTransforms();
  // Output test pictures of two example scenes
  bool TestScenes();
  // Run all tests available
  bool RunAllTests();
};

#endif
