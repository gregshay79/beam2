#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>

struct RGBi {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

RGBi valueToHeatmapColor(double value);
