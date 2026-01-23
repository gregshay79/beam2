#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
//#include <vector>
//#include <iostream>

//#include <windows.h>

// Structure to represent a line segment with color
struct LineSegment {
    float x1, y1, x2, y2;
    float r, g, b; // RGB color

    LineSegment(float x1, float y1, float x2, float y2, float r = 1.0f, float g = 1.0f, float b = 1.0f)
        : x1(x1), y1(y1), x2(x2), y2(y2), r(r), g(g), b(b) {
    }
};

// Window dimensions
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;



class openGLframe
{
public:
    GLuint shaderProgramhandle;
    GLuint VAOhandle, VBOhandle;
    GLFWwindow* window;
    void initLineSegments(std::vector<LineSegment>& _linesegments);
    int setupGL();
    void drawLines(std::vector<LineSegment> _linesegments);
    void waitForCompletion();
    void closeGL();

};

