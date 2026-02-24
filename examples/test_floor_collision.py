"""
Test script to check floor collision forces
"""
import numpy as np
from dem import ParticleSystem

# Create simple system with particle on floor
system = ParticleSystem(
    n_particles=1,
    domain_size=(1.0, 1.0),
    particle_radius=0.1,
    particle_density=2500,
    gravity=9.81,
    k_n=5e4,
    damping=0.8,
    dt=1e-5,
)

# Place particle penetrating into floor
system.positions[0] = [0.5, 0.08]  # Center at y=0.08, bottom at y=-0.02 (penetrating)
system.velocities[0] = [0.0, -0.1]  # Moving downward

print("Testing floor collision:")
print(f"Particle position: {system.positions[0]}")
print(f"Particle velocity: {system.velocities[0]}")
print(f"Particle bottom: y={system.positions[0,1] - system.radius:.3f}m")
print(f"Floor penetration: {system.radius - system.positions[0,1]:.4f}m")
print()

# Compute forces
forces = system.compute_forces()

print(f"Forces on particle: {forces[0]}")
print(f"Gravity force: {system.mass * system.gravity:.2f}N (downward)")
print()

# Expected forces
overlap = system.radius - system.positions[0, 1]
f_n_expected = system.k_n * (overlap ** 1.5)
v_n = system.velocities[0, 1]  # Negative (moving down)
f_d_expected = -system.damping * v_n * np.sqrt(system.k_n * system.mass)

print(f"Expected normal force: {f_n_expected:.2f}N (upward)")
print(f"Expected damping force: {f_d_expected:.2f}N")
print(f"Expected total contact force: {f_n_expected + f_d_expected:.2f}N (upward)")
print(f"Expected net force: {f_n_expected + f_d_expected - system.mass * system.gravity:.2f}N")
print()

# Check if forces are correct
actual_contact_force = forces[0, 1] + system.mass * system.gravity
print(f"Actual contact force: {actual_contact_force:.2f}N")
print(f"Match: {np.isclose(actual_contact_force, f_n_expected + f_d_expected, rtol=0.01)}")
