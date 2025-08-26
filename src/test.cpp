#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <cstring>

#include "yet_another_physics.hpp"  // keeps your particle + physics fns

// ---------------------------- Window & state ----------------------------
static bool isPlaying = true;

static void key_callback(GLFWwindow* w, int key, int sc, int action, int mods) {
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) isPlaying = !isPlaying;
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) glfwSetWindowShouldClose(w, 1);
}

// ---------------------------- Shaders ----------------------------
static const char* kVS = R"(#version 330 core
layout(location=0) in vec2 aQuad;    // per-vertex: unit quad corners in [-1,1]
layout(location=1) in vec2 iCenter;  // per-instance: particle center (world)
layout(location=2) in float iRadius; // per-instance: radius in world units
layout(location=3) in float iVelMag; // per-instance: speed magnitude
layout(location=4) in float iFixed;  // per-instance: 0 or 1

out vec2 vLocal;          // local coords for circle mask
out float vVelMag;
flat out float vFixed;

uniform mat4 uProj;       // orthographic projection

void main(){
    vLocal   = aQuad;           // keep local coords in [-1,1]
    vVelMag  = iVelMag;
    vFixed   = iFixed;

    // Scale quad by radius and translate to particle center
    vec2 world = iCenter + aQuad * iRadius;
    gl_Position = uProj * vec4(world, 0.0, 1.0);
}
)";

static const char* kFS = R"(#version 330 core
in vec2 vLocal;
in float vVelMag;
flat in float vFixed;
out vec4 FragColor;

uniform float uMinVel;
uniform float uMaxVel;

float normVel(float v, float mn, float mx) {
    float denom = max(mx - mn, 1e-6);
    return clamp((v - mn) / denom, 0.0, 1.0);
}

vec3 gradColor(float x){
    float r,g,b;
    if (x < 0.25)      { r=0.0;                g=x*4.0;             b=1.0; }
    else if (x < 0.5)  { r=0.0;                g=1.0;               b=1.0-(x-0.25)*4.0; }
    else if (x < 0.75) { r=(x-0.5)*4.0;        g=1.0;               b=0.0; }
    else               { r=1.0;                g=1.0-(x-0.75)*4.0;  b=0.0; }
    return vec3(r,g,b);
}

void main(){
    // Cut quad into a circle: unit radius in local space
    if (dot(vLocal, vLocal) > 1.0) discard;

    vec3 color = (vFixed > 0.5)
        ? vec3(1.0)                       // fixed particles are white
        : gradColor(normVel(vVelMag, uMinVel, uMaxVel));

    FragColor = vec4(color, 1.0);
}
)";

// ---------------------------- GL helpers ----------------------------
static void checkShader(GLuint s){
    GLint ok=0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if(!ok){ GLint n=0; glGetShaderiv(s, GL_INFO_LOG_LENGTH, &n);
        std::string log(n, '\0'); glGetShaderInfoLog(s, n, &n, log.data());
        std::cerr << "Shader compile error:\n" << log << std::endl; std::exit(1); }
}
static void checkProgram(GLuint p){
    GLint ok=0; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if(!ok){ GLint n=0; glGetProgramiv(p, GL_INFO_LOG_LENGTH, &n);
        std::string log(n, '\0'); glGetProgramInfoLog(p, n, &n, log.data());
        std::cerr << "Program link error:\n" << log << std::endl; std::exit(1); }
}
static GLuint makeProgram(const char* vs, const char* fs){
    GLuint v = glCreateShader(GL_VERTEX_SHADER);   glShaderSource(v,1,&vs,nullptr); glCompileShader(v); checkShader(v);
    GLuint f = glCreateShader(GL_FRAGMENT_SHADER); glShaderSource(f,1,&fs,nullptr); glCompileShader(f); checkShader(f);
    GLuint p = glCreateProgram(); glAttachShader(p,v); glAttachShader(p,f); glLinkProgram(p); checkProgram(p);
    glDeleteShader(v); glDeleteShader(f); return p;
}

// Simple orthographic projection
static void ortho(float l,float r,float b,float t,float n,float f, float outM[16]){
    // Column-major
    outM[0]=2.f/(r-l); outM[4]=0;          outM[8]=0;            outM[12]=-(r+l)/(r-l);
    outM[1]=0;         outM[5]=2.f/(t-b);  outM[9]=0;            outM[13]=-(t+b)/(t-b);
    outM[2]=0;         outM[6]=0;          outM[10]=-2.f/(f-n);  outM[14]=-(f+n)/(f-n);
    outM[3]=0;         outM[7]=0;          outM[11]=0;           outM[15]=1;
}

// ---------------------------- Instance data ----------------------------
struct Instance {
    float cx, cy;     // center
    float radius;     // radius
    float vel;        // speed magnitude
    float fixedFlag;  // 0/1
};

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

int main(){
    // ---- GLFW window/context (Core profile) ----
    if (!glfwInit()){ std::cerr<<"Failed to init GLFW\n"; return -1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    glfwWindowHint(GLFW_OPENGL_PROFILE,GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT,GL_TRUE); // macOS
    GLFWwindow* win = glfwCreateWindow(1200, 800, "SPH Simulation (Shaders)", nullptr, nullptr);
    if(!win){ std::cerr<<"Failed to create window\n"; glfwTerminate(); return -1; }
    glfwMakeContextCurrent(win);
    glfwSetKeyCallback(win, key_callback);
    if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){ std::cerr<<"Failed to load GL\n"; return -1; }
    glClearColor(0.1f,0.1f,0.1f,1.0f);

    // ---- Shader program ----
    GLuint prog = makeProgram(kVS, kFS);
    GLint uProj   = glGetUniformLocation(prog, "uProj");
    GLint uMinVel = glGetUniformLocation(prog, "uMinVel");
    GLint uMaxVel = glGetUniformLocation(prog, "uMaxVel");

    // ---- Geometry: unit quad (two triangles) in local space [-1,1] ----
    float quad[] = {
        -1.f,-1.f,  1.f,-1.f,  1.f, 1.f,
        -1.f,-1.f,  1.f, 1.f, -1.f, 1.f
    };
    GLuint vao=0, vboQuad=0, vboInst=0;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vboQuad);
    glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0);

    // Instance buffer: (cx,cy, radius, vel, fixed)
    glGenBuffers(1, &vboInst);
    glBindBuffer(GL_ARRAY_BUFFER, vboInst);
    glBufferData(GL_ARRAY_BUFFER, 1, nullptr, GL_DYNAMIC_DRAW); // we’ll resize later

    // aCenter (location=1)
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)0);
    glVertexAttribDivisor(1, 1);

    // iRadius (location=2)
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)(2*sizeof(float)));
    glVertexAttribDivisor(2, 1);

    // iVelMag (location=3)
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)(3*sizeof(float)));
    glVertexAttribDivisor(3, 1);

    // iFixed (location=4)
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)(4*sizeof(float)));
    glVertexAttribDivisor(4, 1);

    glBindVertexArray(0);

    // ---- Your simulation setup (unchanged) ----
    double simLastTime = glfwGetTime();
    double fpsLastTime = simLastTime;
    int frameCount = 0;

    write_particles_square_format("../input/particles_2d_5000.txt", 5000, 1.0, 1);

    std::vector<particle> boundary_vec     = create_boundary_particles_2d(1.5, r_e_b);
    std::vector<particle> boundary_circle  = create_boundary_circle_2d(0.5, r_e_b);
    std::vector<particle> boundary_circle_open =
        create_boundary_circle_2d_with_opening(0.5f, 0.04f, 0.2f, M_PI/6.0, M_PI/2.0);

    std::vector<particle> particles(5000);
    readInput("../input/particles_2d_5000.txt", 5000, particles);

    // particles.insert(particles.end(), boundary_circle.begin(), boundary_circle.end());
    particles.insert(particles.end(), boundary_vec.begin(), boundary_vec.end());
    particles.insert(particles.end(), boundary_circle_open.begin(), boundary_circle_open.end());

    std::vector<int> neighbours(particles.size(), -1);
    std::vector<int> particle_neighbours(particles.size(), -1);
    update_neighbors(particles, neighbours, particle_neighbours);

    // ---- Main loop ----
    while(!glfwWindowShouldClose(win)){
        // Time step
        double now = glfwGetTime();
        float dt = 0.001f;
        simLastTime = now;

        if (isPlaying)
            update_physics(dt, particles, neighbours, particle_neighbours);

        // Build instance data + compute min/max velocities
        float min_vel = std::numeric_limits<float>::max();
        float max_vel = std::numeric_limits<float>::lowest();

        std::vector<Instance> inst; inst.reserve(particles.size());
        for (const auto& p : particles){
            float v = std::sqrt(p.vel_x*p.vel_x + p.vel_y*p.vel_y);
            min_vel = std::min(min_vel, v);
            max_vel = std::max(max_vel, v);
        }
        for (const auto& p : particles){
            Instance I;
            I.cx = p.pos_x; I.cy = p.pos_y;
            I.radius = p.radius;
            I.vel = std::sqrt(p.vel_x*p.vel_x + p.vel_y*p.vel_y);
            I.fixedFlag = p.fixed ? 1.0f : 0.0f;
            inst.push_back(I);
        }

        // Upload instance buffer
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        glBufferData(GL_ARRAY_BUFFER, inst.size()*sizeof(Instance), inst.data(), GL_DYNAMIC_DRAW);

        // Set viewport + projection
        int fbw=0, fbh=0; glfwGetFramebufferSize(win,&fbw,&fbh);
        glViewport(0,0,fbw,fbh);
        float aspect = (fbh>0) ? float(fbw)/float(fbh) : 1.0f;

        float P[16];
        if (aspect >= 1.0f)   ortho(-aspect, aspect, -1.f, 1.f, -1.f, 1.f, P);
        else                  ortho(-1.f, 1.f, -1.f/aspect, 1.f/aspect, -1.f, 1.f, P);

        // Draw
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(prog);
        glUniformMatrix4fv(uProj, 1, GL_FALSE, P);
        glUniform1f(uMinVel, min_vel);
        glUniform1f(uMaxVel, max_vel);

        glBindVertexArray(vao);
        glDrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)inst.size());
        glBindVertexArray(0);

        // Title FPS
        frameCount++;
        static double fpsLast = now;
        double elapsed = now - fpsLast;
        if (elapsed >= 0.5){
            double fps = frameCount / elapsed;
            std::string title = "SPH Simulation (Shaders) | t: " + std::to_string(now).substr(0,5)
                              + "s | FPS: " + std::to_string(fps).substr(0,5);
            glfwSetWindowTitle(win, title.c_str());
            frameCount = 0; fpsLast = now;
        }

        glfwSwapBuffers(win);
        glfwPollEvents();
    }

    // Cleanup
    glDeleteBuffers(1, &vboInst);
    glDeleteBuffers(1, &vboQuad);
    glDeleteVertexArrays(1, &vao);
    glDeleteProgram(prog);

    glfwDestroyWindow(win);
    glfwTerminate();
    return 0;
}
