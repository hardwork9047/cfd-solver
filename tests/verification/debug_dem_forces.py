"""
Debug script to check DEM forces
"""
import numpy as np
from dem import ParticleSystem

# Create 4-particle system in smaller domain
system = ParticleSystem(
    n_particles=4,
    domain_size=(0.6, 0.8),
    particle_radius=0.1,
    particle_density=2500,
    gravity=9.81,
    k_n=1.33e4,
    damping=0.8,
    dt=1e-5,
)

# Place particles without initial overlap - they will contact when falling
system.positions[0] = [0.15, 0.12]  # Bottom left
system.positions[1] = [0.35, 0.12]  # Bottom right (closer)
system.positions[2] = [0.15, 0.25]  # Top left - just above particle 0 (will contact soon)
system.positions[3] = [0.35, 0.25]  # Top right - just above particle 1 (will contact soon)
system.velocities[:] = 0.0

print("Initial state:")
for i in range(system.n_particles):
    print(f"Particle {i}: pos={system.positions[i]}, vel={system.velocities[i]}")
print()

# Check distances between all pairs
print("Initial distances:")
for i in range(system.n_particles):
    for j in range(i+1, system.n_particles):
        dist = np.linalg.norm(system.positions[j] - system.positions[i])
        overlap = 2*system.radius - dist
        print(f"Particles {i}-{j}: distance={dist:.4f}m, overlap={overlap:.4f}m")
print()

# Compute forces
forces = system.compute_forces()

print("Forces:")
for i in range(system.n_particles):
    print(f"Particle {i}: {forces[i]}")
print()

print(f"Particle mass: {system.mass:.2f}kg")
print(f"Gravity force per particle: {system.mass * system.gravity:.2f}N")
print()

# Calculate expected force for first overlap
overlap = 2*system.radius - np.linalg.norm(system.positions[1] - system.positions[0])
f_n_expected = system.k_n * (overlap ** 1.5) if overlap > 0 else 0
gravity_force = system.mass * system.gravity

print(f"Example normal force (particles 0-1): {f_n_expected:.2f}N")
print(f"Gravity force: {gravity_force:.2f}N")
if f_n_expected > 0:
    print(f"Ratio (normal/gravity): {f_n_expected/gravity_force:.4f}")
print()

# Time evolution
print("Running simulation for 2000 steps...")
for step in range(2000):
    system.step()
    
    if step % 400 == 0 or step == 1999:
        print(f"\nStep {step} (t={step*system.dt:.4f}s):")
        for i in range(system.n_particles):
            print(f"  Particle {i}: pos={system.positions[i]}, vel={system.velocities[i]}")
        
        # Check overlaps
        overlap_count = 0
        for i in range(system.n_particles):
            for j in range(i+1, system.n_particles):
                distance = np.linalg.norm(system.positions[j] - system.positions[i])
                overlap = 2*system.radius - distance
                if overlap > 0:
                    overlap_count += 1
                    print(f"  Particles {i}-{j}: dist={distance:.4f}m, overlap={overlap:.4f}m")
        if overlap_count == 0:
            print("  No overlaps")
