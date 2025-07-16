import numpy as np
import random
import os

# === Setup: run from the script's own directory ===
path_Script = os.path.dirname(os.path.realpath(__file__))
os.chdir(path_Script)

# === VTK Writer Function ===
def write_vtk_file(filename, block):
    with open(filename, 'w') as vtk_file:
        vtk_file.write('# vtk DataFile Version 4.0\n')
        vtk_file.write('Particle data\n')
        vtk_file.write('ASCII\n')
        vtk_file.write('DATASET POLYDATA\n')

        num_particles = len(block)
        vtk_file.write(f'POINTS {num_particles} float\n')
        for p in block:
            vtk_file.write(f"{p['x']} {p['y']} {p['z']}\n")

        vtk_file.write(f'VERTICES {num_particles} {2 * num_particles}\n')
        for i in range(num_particles):
            vtk_file.write(f'1 {i}\n')

        vtk_file.write(f'POINT_DATA {num_particles}\n')

        vtk_file.write('SCALARS mass double\nLOOKUP_TABLE default\n')
        for p in block:
            vtk_file.write(f"{p['mass']}\n")

        vtk_file.write('VECTORS velocity double\n')
        for p in block:
            vtk_file.write(f"{p['vx']} {p['vy']} {p['vz']}\n")

        vtk_file.write('SCALARS radius double\nLOOKUP_TABLE default\n')
        for p in block:
            vtk_file.write(f"{p['radius']}\n")

        vtk_file.write('SCALARS fixed int\nLOOKUP_TABLE default\n')
        for p in block:
            vtk_file.write(f"{p['fixed']}\n")

# === Wall particle generator ===
def generate_wall_particles():
    """
    Generate wall particles forming an open-top box.
    All particles are fixed (fixed = 1) and stationary.
    """
    particles = []
    dx = dy = dz = 0.0075

    x_min, x_max = 1.125, 1.875
    y_min, y_max = 1.5, 2.25
    z_min, z_max = 1.125, 1.875

    x_start, y_start, z_start = 1.10, 1.45, 1.10
    Nx = int((1.90 - x_start) / dx) + 1
    Ny = int((2.30 - y_start) / dy) + 1
    Nz = int((1.90 - z_start) / dz) + 1

    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_start + dx * i
                y = y_start + dy * j
                z = z_start + dz * k

                if x < x_min or x > x_max or y < y_min or y > y_max or z < z_min or z > z_max:
                    if y >= y_max:
                        continue  # Open-top box
                    particle = {
                        'x': x, 'y': y, 'z': z,
                        'vx': 0, 'vy': 0, 'vz': 0,
                        'mass': 0.012, 'radius': 0.01,
                        'fixed': 1 
                    }
                    particles.append(particle)

    return particles

# === Fluid particle generator ===
def generate_particles():
    particles = []
    keep_ratio = 0.3  # Keep 10%
    radius = 0.3      # Tiny fluid blob

    for i in range(51):
        for j in range(51):
            for k in range(51):
                x = 1.0 + 0.02 * i
                y = 1.99 + 0.02 * j
                z = 1.0 + 0.02 * k

                if np.linalg.norm(np.array([x, y, z]) - np.array([1.5, 2.49, 1.5])) <= radius:
                    if random.random() <= keep_ratio:
                        particle = {
                            'x': x,
                            'y': y,
                            'z': z,
                            'vx': 0,
                            'vy': 0,
                            'vz': 0,
                            'mass': 0.012,
                            'radius': 0.01,
                            'fixed': 0
                        }
                        particles.append(particle)
    return particles

def count_particles(particles):
    total = len(particles)
    fluid = sum(1 for p in particles if p['fixed'] == 0)
    wall = total - fluid
    print(f"Total particles: {total}")
    print(f"  Fluid particles: {fluid}")
    print(f"  Wall particles:  {wall}")


# === Main execution block ===
if __name__ == "__main__":
    fluid_particles = generate_particles()
    wall_particles = generate_wall_particles()
    all_particles = fluid_particles + wall_particles
    count_particles(all_particles)
    write_vtk_file("particle_scene.vtk", all_particles)