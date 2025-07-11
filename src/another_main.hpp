#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "physics_debug.hpp"
#include "utils/dec.hpp"
#include "utils/particle.hpp"
#include "utils/neighbour.hpp"

const int NUM_SEGMENTS = 100;
const float RADIUS = 0.05f;

bool isPlaying = true;
// bool applyClickForce = false;
// float click_x = 0.0f, click_y = 0.0f;
bool isClickActive = false; // true while right mouse button is held
float click_x = 0.0f, click_y = 0.0f;

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if (action == GLFW_PRESS)
        {
            isClickActive = true;
        }
        else if (action == GLFW_RELEASE)
        {
            isClickActive = false;
        }
    }
}

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
    {
        isPlaying = !isPlaying;
    }
}

void apply_click_repulsion(std::vector<particle> &particles, float cx, float cy, float dt)
{
    const float strength = 50.0f;
    const float effective_radius = 0.4f;

    for (auto &p : particles)
    {
        if (p.is_boundary)
            continue;
        float dx = p.x - cx;
        float dy = p.y - cy;
        float r2 = dx * dx + dy * dy;

        if (r2 > effective_radius * effective_radius)
            continue;

        float r = std::sqrt(r2);
        if (r < 1e-5f)
            r = 1e-5f;

        float force_mag = strength * (1.0f - r / effective_radius);
        float fx = (dx / r) * force_mag;
        float fy = (dy / r) * force_mag;

        p.velocity_x += fx * dt;
        p.velocity_y += fy * dt;
    }
}

void drawCircle(float cx, float cy, float r, int segments)
{
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(cx, cy); // center of circle
    for (int i = 0; i <= segments; ++i)
    {
        float theta = 2.0f * M_PI * float(i) / float(segments);
        float x = r * cosf(theta);
        float y = r * sinf(theta);
        glVertex2f(x + cx, y + cy);
    }
    glEnd();
}

void update_physics(std::vector<particle>&particles,std::vector<std::vector<int>>&buckets,int nBuckets,int boundary,sph_float dt)
{
    // std::cout<<"the dt is :"<<dt<<std::endl;
    // Update physics
    update_buckets(particles, buckets, boundary, nBuckets);
    // calculate_density(particles);
    calculate_density_buckets(particles, buckets, nBuckets);
    calculate_pressure(particles);
    //calculate_force(particles, dt);
    calculate_force_buckets(particles, buckets, nBuckets, dt);
    //     if (isClickActive)
    // {
    //     double xpos, ypos;
    //     int width, height;
    //     glfwGetCursorPos(glfwGetCurrentContext(), &xpos, &ypos);
    //     glfwGetFramebufferSize(glfwGetCurrentContext(), &width, &height);

    //     float x_ndc = 2.0f * xpos / width - 1.0f;
    //     float y_ndc = 1.0f - 2.0f * ypos / height;

    //     float aspect = static_cast<float>(width) / height;
    //     float cx = (aspect >= 1.0f) ? x_ndc * aspect : x_ndc;
    //     float cy = (aspect >= 1.0f) ? y_ndc : y_ndc / aspect;

    //     apply_click_repulsion(particles, cx, cy, dt);
    // }

    // apply_gravity(particles, -9.8, dt);
    update_positions(particles, dt);
    // apply_boundary_conditions(particles, 2.0f, 2.0f);
    //  std::cout<<particles.size()<<" - is the particles vector size"<<std::endl;
}

// Gradient: blue → cyan → green → yellow → red
void velocityToColor(float v, float &r, float &g, float &b)
{
    v = std::clamp(v, 0.0f, 1.0f); // Ensure within bounds

    if (v < 0.25f)
    {
        // Blue → Cyan
        r = 0.0f;
        g = v * 4.0f;
        b = 1.0f;
    }
    else if (v < 0.5f)
    {
        // Cyan → Green
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - (v - 0.25f) * 4.0f;
    }
    else if (v < 0.75f)
    {
        // Green → Yellow
        r = (v - 0.5f) * 4.0f;
        g = 1.0f;
        b = 0.0f;
    }
    else
    {
        // Yellow → Red
        r = 1.0f;
        g = 1.0f - (v - 0.75f) * 4.0f;
        b = 0.0f;
    }
}

void drawScene(std::vector<particle>&particles,GLFWwindow *window)
{
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    float aspect = float(width) / float(height);
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (aspect >= 1.0f)
        glOrtho(-aspect, aspect, -1, 1, -1, 1);
    else
        glOrtho(-1, 1, -1 / aspect, 1 / aspect, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Render particles
    glClear(GL_COLOR_BUFFER_BIT);
    // 1. Compute min and max velocity magnitudes
    float min_vel = std::numeric_limits<float>::max();
    float max_vel = std::numeric_limits<float>::lowest();

    for (const auto &p : particles)
    {
        float v = std::sqrt(p.velocity_x * p.velocity_x + p.velocity_y * p.velocity_y);
        min_vel = std::min(min_vel, v);
        max_vel = std::max(max_vel, v);
    }

    // 2. Draw particles with color gradient
    for (const auto &p : particles)
    {
        if (p.is_boundary)
        {
            glColor3f(1.0f, 1.0f, 1.0f); // White color
        }
        else
        {
            float v = std::sqrt(p.velocity_x * p.velocity_x + p.velocity_y * p.velocity_y);
            float norm_v = (v - min_vel) / (max_vel - min_vel + 1e-5f);

            float r, g, b;
            velocityToColor(norm_v, r, g, b);
            glColor3f(r, g, b);
        }

        drawCircle(p.x, p.y, p.radius, NUM_SEGMENTS);
    }
}

int main()
{
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    GLFWwindow *window = glfwCreateWindow(1200, 800, "SPH Simulation", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        std::cerr << "Failed to create window\n";
        return -1;
    }

    glfwMakeContextCurrent(window);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f); // background color

    // Orthographic 2D projection (updated dynamically)
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);

    // Initialize particles in a square

    double simLastTime = glfwGetTime();
    double fpsLastTime = simLastTime;
    int frameCount = 0;

    // initialization
    std::vector<particle> particles = create_all_particles();
    int boundary = 3;
    int nBuckets = ceil(boundary / r_e);
    std::vector<std::vector<int>> buckets(nBuckets * nBuckets);
    std::vector<int> neighbour_list(particles.size());
    std::vector<int> start_neighbour(particles.size(), -1);
    //init

    update_buckets(particles, buckets, boundary, nBuckets);

    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    while (!glfwWindowShouldClose(window))
    {
        double currentTime = glfwGetTime();
        float dt = static_cast<float>(currentTime - simLastTime);
        simLastTime = currentTime;

        if (isPlaying)
            update_physics(particles,buckets,boundary,nBuckets,dt);
        // Update viewport and aspect ratio
        drawScene(particles,window);
        glfwSwapBuffers(window);
        glfwPollEvents();

        // FPS and time display
        frameCount++;
        double elapsedTime = currentTime - fpsLastTime;
        if (elapsedTime >= 0.5)
        {
            double fps = frameCount / elapsedTime;
            std::string title = "SPH Simulation | Time: " + std::to_string(currentTime).substr(0, 5) +
                                "s | FPS: " + std::to_string(fps).substr(0, 5);
            glfwSetWindowTitle(window, title.c_str());

            frameCount = 0;
            fpsLastTime = currentTime;
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
