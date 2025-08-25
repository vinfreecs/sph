#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "yet_another_physics.hpp"

const int NUM_SEGMENTS = 100;
const float RADIUS = 0.05f;

bool isPlaying = true;
bool isClickActive = false;     // true while right mouse button is held
float click_x = 0.0f, click_y = 0.0f;

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
    {
        isPlaying = !isPlaying;
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

void update_physics(float dt,std::vector<particle>&particles,std::vector<int>& neighbours,std::vector<int>& particle_neighbours)
{
    // std::cout<<"the dt is :"<<dt<<std::endl;
    // Update physics
    update_neighbors(particles, neighbours, particle_neighbours);
    calculation_density_pressure_hash(particles, neighbours, particle_neighbours);
    reset_forces(particles);
    calculate_force_hash(particles, neighbours, particle_neighbours);
    update_velocity(particles, dt);
    update_positions(particles, dt);
    // apply_boundary_conditions(particles, -1.0f, 1.0f, -1.0f, 1.0f);
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

void drawScene(GLFWwindow *window,std::vector<particle>& particles)
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
        float v = std::sqrt(p.vel_x * p.vel_x + p.vel_y * p.vel_y);
        min_vel = std::min(min_vel, v);
        max_vel = std::max(max_vel, v);
    }

    // 2. Draw particles with color gradient
    for (const auto &p : particles)
    {
        if(!p.fixed){
            float v = std::sqrt(p.vel_x * p.vel_x + p.vel_y * p.vel_y);
            float norm_v = (v - min_vel) / (max_vel - min_vel + 1e-5f);
    
            float r, g, b;
            velocityToColor(norm_v, r, g, b);
    
            glColor3f(r, g, b);
            drawCircle(p.pos_x, p.pos_y, p.radius, NUM_SEGMENTS);

        }
        else{
            glColor3f(255, 255, 255);
            drawCircle(p.pos_x, p.pos_y, p.radius, NUM_SEGMENTS);
        }
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
    write_particles_square_format("../input/particles_2d_1000.txt", 1000, 1.0,1);

    std::vector<particle> boundary_vec = create_boundary_particles_2d(1.5,r_e_b);
    std::vector<particle> boundary_circle = create_boundary_circle_2d(0.5,r_e_b);
    std::vector<particle> boundary_circle_open =create_boundary_circle_2d_with_opening(0.5f, 0.04f, 0.2f, M_PI / 6.0, M_PI / 2.0);




    std::vector<particle> particles(1000);
    readInput("../input/particles_2d_1000.txt",1000,particles);

    // particles.insert(particles.end(),boundary_circle.begin(),boundary_circle.end());
    particles.insert(particles.end(),boundary_vec.begin(),boundary_vec.end());
    particles.insert(particles.end(),boundary_circle_open.begin(),boundary_circle_open.end());

    std::vector<int> neighbours(particles.size(), -1);
    std::vector<int> particle_neighbours(particles.size(),-1);

    update_neighbors(particles, neighbours, particle_neighbours);

    glfwSetKeyCallback(window, key_callback);
    while (!glfwWindowShouldClose(window))
    {
        double currentTime = glfwGetTime();
        float dt = static_cast<float>(currentTime - simLastTime);
        dt = 0.001;
        simLastTime = currentTime;

        if (isPlaying)
            update_physics(dt,particles, neighbours,particle_neighbours);
        // Update viewport and aspect ratio
        drawScene(window,particles);
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
