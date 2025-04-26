#include <iostream>
#include <cmath>
#include <algorithm>

struct RGB {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

RGB valueToHeatmapColor(double value);
