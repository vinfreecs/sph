// 2D OpenGL viewer: runs the shared CPU solver on the "cup" scene and draws particles as
// instanced circles coloured by speed (blue = slow, red = fast). Space = pause, Esc = quit.

#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "sph/app.hpp"
#include "sph/cpu_solver.hpp"

static bool g_playing = true;

static void key_callback(GLFWwindow* w, int key, int, int action, int) {
    if (action != GLFW_PRESS) return;
    if (key == GLFW_KEY_SPACE) g_playing = !g_playing;
    if (key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(w, 1);
}

static const char* kVS = R"(#version 330 core
layout(location=0) in vec2 aQuad;     // unit quad corner in [-1,1]
layout(location=1) in vec2 iCenter;   // per instance: particle centre
layout(location=2) in float iRadius;
layout(location=3) in float iSpeed;
layout(location=4) in float iFixed;
out vec2 vLocal;
out float vSpeed;
flat out float vFixed;
uniform mat4 uProj;
void main() {
    vLocal = aQuad;
    vSpeed = iSpeed;
    vFixed = iFixed;
    gl_Position = uProj * vec4(iCenter + aQuad * iRadius, 0.0, 1.0);
}
)";

static const char* kFS = R"(#version 330 core
in vec2 vLocal;
in float vSpeed;
flat in float vFixed;
out vec4 FragColor;
uniform float uMinSpeed;
uniform float uMaxSpeed;
vec3 ramp(float x) {               // blue -> cyan -> green -> yellow -> red
    if (x < 0.25) return vec3(0.0, x * 4.0, 1.0);
    if (x < 0.50) return vec3(0.0, 1.0, 1.0 - (x - 0.25) * 4.0);
    if (x < 0.75) return vec3((x - 0.5) * 4.0, 1.0, 0.0);
    return vec3(1.0, 1.0 - (x - 0.75) * 4.0, 0.0);
}
void main() {
    if (dot(vLocal, vLocal) > 1.0) discard;
    float t = clamp((vSpeed - uMinSpeed) / max(uMaxSpeed - uMinSpeed, 1e-6), 0.0, 1.0);
    FragColor = vec4(vFixed > 0.5 ? vec3(1.0) : ramp(t), 1.0);
}
)";

static GLuint compile(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetShaderInfoLog(s, sizeof log, nullptr, log);
        std::cerr << "shader compile error:\n" << log << "\n";
        std::exit(1);
    }
    return s;
}

static GLuint make_program() {
    GLuint v = compile(GL_VERTEX_SHADER, kVS), f = compile(GL_FRAGMENT_SHADER, kFS);
    GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetProgramInfoLog(p, sizeof log, nullptr, log);
        std::cerr << "program link error:\n" << log << "\n";
        std::exit(1);
    }
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

static void ortho(float l, float r, float b, float t, float m[16]) {  // column-major, z in [-1,1]
    std::fill(m, m + 16, 0.0f);
    m[0] = 2 / (r - l);
    m[5] = 2 / (t - b);
    m[10] = -1;
    m[12] = -(r + l) / (r - l);
    m[13] = -(t + b) / (t - b);
    m[15] = 1;
}

struct Instance {
    float cx, cy, radius, speed, fixed;
};

int main(int argc, char** argv) {
    int fluid = argc > 1 ? std::atoi(argv[1]) : 1000;
    sph::Scene scene;
    if (!sph::make_scene("cup", fluid, scene)) return 1;
    sph::print_scene_summary(scene);
    sph::CpuSolver<2> solver;

    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);  // required on macOS
    GLFWwindow* win = glfwCreateWindow(1200, 800, "SPH (OpenGL)", nullptr, nullptr);
    if (!win) {
        std::cerr << "failed to create window\n";
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(win);
    glfwSetKeyCallback(win, key_callback);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "failed to load OpenGL\n";
        return 1;
    }
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

    GLuint prog = make_program();
    GLint uProj = glGetUniformLocation(prog, "uProj");
    GLint uMin = glGetUniformLocation(prog, "uMinSpeed");
    GLint uMax = glGetUniformLocation(prog, "uMaxSpeed");

    const float quad[] = {-1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1};
    GLuint vao, vbo_quad, vbo_inst;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo_quad);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_quad);
    glBufferData(GL_ARRAY_BUFFER, sizeof quad, quad, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);

    glGenBuffers(1, &vbo_inst);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_inst);
    const int sizes[] = {2, 1, 1, 1};
    size_t offset = 0;
    for (int a = 0; a < 4; ++a) {
        glEnableVertexAttribArray(a + 1);
        glVertexAttribPointer(a + 1, sizes[a], GL_FLOAT, GL_FALSE, sizeof(Instance), (void*)offset);
        glVertexAttribDivisor(a + 1, 1);
        offset += sizes[a] * sizeof(float);
    }
    glBindVertexArray(0);

    std::vector<Instance> inst;
    double fps_t0 = glfwGetTime();
    int frames = 0;
    while (!glfwWindowShouldClose(win)) {
        if (g_playing) solver.step(scene.particles, scene.params, scene.dt);

        float vmin = std::numeric_limits<float>::max(), vmax = 0;
        inst.clear();
        for (const auto& p : scene.particles) {
            float v = std::sqrt(p.vel_x * p.vel_x + p.vel_y * p.vel_y);
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
            inst.push_back({p.pos_x, p.pos_y, p.radius, v, p.fixed ? 1.0f : 0.0f});
        }
        glBindBuffer(GL_ARRAY_BUFFER, vbo_inst);
        glBufferData(GL_ARRAY_BUFFER, inst.size() * sizeof(Instance), inst.data(), GL_DYNAMIC_DRAW);

        int w, h;
        glfwGetFramebufferSize(win, &w, &h);
        glViewport(0, 0, w, h);
        float aspect = h > 0 ? float(w) / h : 1.0f, P[16];
        if (aspect >= 1) ortho(-aspect, aspect, -1, 1, P);
        else ortho(-1, 1, -1 / aspect, 1 / aspect, P);

        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(prog);
        glUniformMatrix4fv(uProj, 1, GL_FALSE, P);
        glUniform1f(uMin, vmin);
        glUniform1f(uMax, vmax);
        glBindVertexArray(vao);
        glDrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)inst.size());
        glfwSwapBuffers(win);
        glfwPollEvents();

        ++frames;
        double now = glfwGetTime();
        if (now - fps_t0 >= 0.5) {
            std::string title = "SPH (OpenGL) | " + std::to_string(int(frames / (now - fps_t0))) + " FPS" +
                                (g_playing ? "" : " | paused");
            glfwSetWindowTitle(win, title.c_str());
            frames = 0;
            fps_t0 = now;
        }
    }
    glDeleteBuffers(1, &vbo_inst);
    glDeleteBuffers(1, &vbo_quad);
    glDeleteVertexArrays(1, &vao);
    glDeleteProgram(prog);
    glfwDestroyWindow(win);
    glfwTerminate();
}
