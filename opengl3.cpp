#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <iostream>

#include <windows.h>
#include "opengl3.h"

// Vector to store line segments
//std::vector<LineSegment> lineSegments;


// Error callback function
void errorCallback(int error, const char* description) {
    std::cerr << "GLFW Error " << error << ": " << description << std::endl;
}

// Key callback function
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
}



    // Initialize line segments for demonstration
    void mygraphics::initLineSegments(std::vector<LineSegment> &_lineSegments) {
        // Add some sample line segments with different colors
        _lineSegments.emplace_back(-0.8f, -0.8f, 0.8f, 0.8f, 1.0f, 0.0f, 0.0f);  // Red diagonal
        _lineSegments.emplace_back(0.8f, -0.8f, -0.8f, 0.8f, 0.0f, 1.0f, 0.0f);  // Green diagonal
        _lineSegments.emplace_back(-0.5f, 0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f);   // Blue horizontal
        _lineSegments.emplace_back(0.0f, -0.5f, 0.0f, 0.5f, 1.0f, 1.0f, 0.0f);   // Yellow vertical

        // Add a square
        float s = 0.3f;
        _lineSegments.emplace_back(-s, -s, s, -s, 1.0f, 0.5f, 0.0f);  // Bottom
        _lineSegments.emplace_back(s, -s, s, s, 1.0f, 0.5f, 0.0f);    // Right
        _lineSegments.emplace_back(s, s, -s, s, 1.0f, 0.5f, 0.0f);    // Top
        _lineSegments.emplace_back(-s, s, -s, -s, 1.0f, 0.5f, 0.0f);  // Left
    }

    int mygraphics::setupGL() {
        // Initialize GLFW
        if (!glfwInit()) {
            std::cerr << "Failed to initialize GLFW" << std::endl;
            return -1;
        }

        // Set error callback
        glfwSetErrorCallback(errorCallback);

        // Set GLFW window hints
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        // Create a windowed mode window and its OpenGL context
        window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "OpenGL Line Segment Plotter", NULL, NULL);
        if (!window) {
            std::cerr << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            return -1;
        }

        // Make the window's context current
        glfwMakeContextCurrent(window);

        // Set key callback
        glfwSetKeyCallback(window, keyCallback);

        // Initialize GLEW
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) {
            std::cerr << "Failed to initialize GLEW" << std::endl;
            glfwDestroyWindow(window);
            glfwTerminate();
            return -1;
        }

        // Create and compile vertex shader
        const char* vertexShaderSource =
            "#version 330 core\n"
            "layout (location = 0) in vec2 position;\n"
            "layout (location = 1) in vec3 color;\n"
            "out vec3 vertexColor;\n"
            "void main() {\n"
            "    gl_Position = vec4(position.x, position.y, 0.0, 1.0);\n"
            "    vertexColor = color;\n"
            "}\0";

        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glCompileShader(vertexShader);

        // Check for vertex shader compile errors
        int success;
        char infoLog[512];
        glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
            std::cerr << "Vertex shader compilation failed: " << infoLog << std::endl;
        }

        // Create and compile fragment shader
        const char* fragmentShaderSource =
            "#version 330 core\n"
            "in vec3 vertexColor;\n"
            "out vec4 fragColor;\n"
            "void main() {\n"
            "    fragColor = vec4(vertexColor, 1.0);\n"
            "}\0";

        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(fragmentShader);

        // Check for fragment shader compile errors
        glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
            std::cerr << "Fragment shader compilation failed: " << infoLog << std::endl;
        }

        // Link shaders
        shaderProgramhandle = glCreateProgram();
        glAttachShader(shaderProgramhandle, vertexShader);
        glAttachShader(shaderProgramhandle, fragmentShader);
        glLinkProgram(shaderProgramhandle);

        // Check for linking errors
        glGetProgramiv(shaderProgramhandle, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(shaderProgramhandle, 512, NULL, infoLog);
            std::cerr << "Shader program linking failed: " << infoLog << std::endl;
        }

        // Delete the shaders as they're linked into our program now and no longer necessary
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        // Initialize line segments
        //initLineSegments();



        // Set up vertex buffer object (VBOhandle) and vertex array object (VAOhandle)
        //GLuint VAOhandle, VBOhandle;
        glGenVertexArrays(1, &VAOhandle);
        glGenBuffers(1, &VBOhandle);

        // Bind the VAOhandle first, then bind and set VBOhandle, and then configure vertex attributes
        //glBindVertexArray(VAOhandle);

        //glBindBuffer(GL_ARRAY_BUFFER, VBOhandle);
        //glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

        //// Position attribute (x, y)
        //glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        //glEnableVertexAttribArray(0);

        //// Color attribute (r, g, b)
        //glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(2 * sizeof(float)));
        //glEnableVertexAttribArray(1);

        //// Unbind VBOhandle and VAOhandle
        //glBindBuffer(GL_ARRAY_BUFFER, 0);
        //glBindVertexArray(0);

        // Main loop
        //while (!glfwWindowShouldClose(window)) {
            // Clear the screen
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Use shader program
        //glUseProgram(shaderProgramhandle);

        // Bind VAOhandle
        //glBindVertexArray(VAOhandle);

        // Select shader program to use
        glUseProgram(shaderProgramhandle);

        // Bind the VAOhandle first, then bind and set VBOhandle, and then configure vertex attributes
        glBindVertexArray(VAOhandle);
        glBindBuffer(GL_ARRAY_BUFFER, VBOhandle);

        // Position attribute (x, y)
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        // Color attribute (r, g, b)
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);

        return 0;
    }

    void mygraphics::draw(std::vector<LineSegment> _lineSegments)
    {

        // Prepare vertex data for all line segments
        std::vector<float> vertices;

        // Each vertex needs 5 floats: x, y, r, g, b
        for (const auto& line : _lineSegments) {
            // First point of the line
            vertices.push_back(line.x1);    // x position
            vertices.push_back(line.y1);    // y position
            vertices.push_back(line.r);     // red
            vertices.push_back(line.g);     // green
            vertices.push_back(line.b);     // blue

            // Second point of the line
            vertices.push_back(line.x2);    // x position
            vertices.push_back(line.y2);    // y position
            vertices.push_back(line.r);     // red
            vertices.push_back(line.g);     // green
            vertices.push_back(line.b);     // blue
        }

        // Clear the screen
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        //// Use shader program
        //glUseProgram(shaderProgramhandle);
        //// Bind the VAOhandle first, then bind and set VBOhandle, and then configure vertex attributes
        //glBindVertexArray(VAOhandle);
        //glBindBuffer(GL_ARRAY_BUFFER, VBOhandle);

        //Shader, VAO and VBO are already held bound.
        // Vertex attributes are already set.
        // just configure buffer data with vertex data here.
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

        //// Position attribute (x, y)
        //glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        //glEnableVertexAttribArray(0);

        //// Color attribute (r, g, b)
        //glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(2 * sizeof(float)));
        //glEnableVertexAttribArray(1);

        // Draw all line segments in a single call
        glDrawArrays(GL_LINES, 0, (GLsizei)(_lineSegments.size() * 2));

        //// Unbind VBOhandle and VAOhandle
        //glBindBuffer(GL_ARRAY_BUFFER, 0);
        //glBindVertexArray(0);

        // Swap front and back buffers
        glfwSwapBuffers(window);
    }

    // Poll for and process events
    void mygraphics::waitForCompletion() {
        while (!glfwWindowShouldClose(window)) {

            glfwPollEvents();
        }
    }

    void mygraphics::closeGL()
    {
        // Unbind VBOhandle and VAOhandle
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        // Clean up resources
        glDeleteVertexArrays(1, &VAOhandle);
        glDeleteBuffers(1, &VBOhandle);
        glDeleteProgram(shaderProgramhandle);

        // Terminate GLFW
        glfwDestroyWindow(window);
        glfwTerminate();
    }


//int main2()
//{
//    mygraphics graphics;
//
//    graphics.setupGL();
//
//    // Initialize line segments
//    graphics.initLineSegments();
//
//    //graphics.runGL();
//
//    //graphics.prepareDrawing();
//    graphics.draw(lineSegments);
//
//    for (int i = 0; i < 1000; i++) {
//        //Sleep(100);
//        for (auto& line : lineSegments) {
//            line.x1 *= .99;
//            line.x2 *= .99;
//            line.y1 *= .99;
//            line.y2 *= .99;
//        }
//        graphics.draw(lineSegments);
//
//        glfwPollEvents();
//        if (glfwWindowShouldClose(graphics.window))
//            break;
//    }
//
//    graphics.waitForCompletion();
//
//    graphics.closeGL();
//
//    return 0;
//}