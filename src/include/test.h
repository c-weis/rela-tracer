// Copyright 2023 Christoph Weis
#ifndef TEST_H_INCLUDED
#define TEST_H_INCLUDED

#include "tracer.h"

class Tester {
 public:
    bool TestRandomness();
    bool TestColors();
    bool TestLineLorentzTransforms();
    bool TestColorLorentzTransforms();
    bool TestScenes();
    bool RunAllTests();
};

#endif