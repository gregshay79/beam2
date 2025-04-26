#include <iostream>
#include <cmath>
#include <algorithm>

struct RGB {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

RGB valueToHeatmapColor(double value) {
    // Ensure value is in the range [0.0, 1.0]
    value = std::clamp(value, 0.0, 1.0);

    RGB color;

    // Classic heat map: blue -> cyan -> green -> yellow -> red
    if (value < 0.25) {
        // Blue to Cyan (0.0 - 0.25)
        double normalizedValue = value / 0.25;
        color.r = 0;
        color.g = static_cast<unsigned char>(255 * normalizedValue);
        color.b = 255;
    }
    else if (value < 0.5) {
        // Cyan to Green (0.25 - 0.5)
        double normalizedValue = (value - 0.25) / 0.25;
        color.r = 0;
        color.g = 255;
        color.b = static_cast<unsigned char>(255 * (1.0 - normalizedValue));
    }
    else if (value < 0.75) {
        // Green to Yellow (0.5 - 0.75)
        double normalizedValue = (value - 0.5) / 0.25;
        color.r = static_cast<unsigned char>(255 * normalizedValue);
        color.g = 255;
        color.b = 0;
    }
    else {
        // Yellow to Red (0.75 - 1.0)
        double normalizedValue = (value - 0.75) / 0.25;
        color.r = 255;
        color.g = static_cast<unsigned char>(255 * (1.0 - normalizedValue));
        color.b = 0;
    }

    return color;
}

// Example usage
void testHeatmap() {
    std::cout << "Heatmap Color Examples:" << std::endl;
    for (double value = 0.0; value <= 1.0; value += 0.1) {
        RGB color = valueToHeatmapColor(value);
        std::cout << "Value: " << value << " -> RGB: ("
            << static_cast<int>(color.r) << ", "
            << static_cast<int>(color.g) << ", "
            << static_cast<int>(color.b) << ")" << std::endl;
    }
}