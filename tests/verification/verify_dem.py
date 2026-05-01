"""
Verification script for DEM particle system.

This script demonstrates and validates the DEM implementation.
"""

import matplotlib.pyplot as plt
import numpy as np

from dem import ParticleSystem


def verify_single_particle_fall():
    """Verify single particle falling under gravity."""
    print("=" * 60)
    print("DEM Verification - Single Particle Free Fall")
    print("=" * 60)

    system = ParticleSystem(
        n_particles=1,
        domain_size=(5.0, 10.0),
        particle_radius=0.2,
        particle_density=2500,
        gravity=9.81,
        k_n=1e6,
        damping=0.1,
        dt=1e-4,
    )

    # Set initial position at top center
    system.positions[0] = [2.5, 9.0]
    system.velocities[0] = [0.0, 0.0]

    print(f"\nInitial position: {system.positions[0]}")
    print(f"Initial velocity: {system.velocities[0]}")
    print(f"Particle mass: {system.mass:.3e}kg")

    # Record trajectory
    times = []
    heights = []
    velocities_y = []

    t_end = 1.5
    n_steps = int(t_end / system.dt)

    for _ in range(n_steps):
        times.append(system.time)
        heights.append(system.positions[0, 1])
        velocities_y.append(system.velocities[0, 1])
        system.step()

    print(f"\nFinal position: {system.positions[0]}")
    print(f"Final velocity: {system.velocities[0]}")

    # Plot trajectory
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(times, heights, "b-", linewidth=2, label="Simulated")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Height [m]")
    ax1.set_title("Particle Height vs Time")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(times, velocities_y, "r-", linewidth=2, label="Simulated")
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Vertical Velocity [m/s]")
    ax2.set_title("Particle Velocity vs Time")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("dem_single_particle_fall.png", dpi=150, bbox_inches="tight")
    print("\nPlot saved as: dem_single_particle_fall.png")

    plt.show()


def verify_particle_collision():
    """Verify two-particle collision."""
    print("\n" + "=" * 60)
    print("DEM Verification - Two Particle Collision")
    print("=" * 60)

    system = ParticleSystem(
        n_particles=2,
        domain_size=(5.0, 5.0),
        particle_radius=0.3,
        particle_density=2500,
        gravity=0.0,  # No gravity for collision test
        k_n=1e5,
        damping=0.2,
        dt=1e-4,
    )

    # Set initial positions (approaching each other)
    system.positions[0] = [1.5, 2.5]
    system.positions[1] = [3.5, 2.5]
    system.velocities[0] = [2.0, 0.0]
    system.velocities[1] = [-2.0, 0.0]

    print(f"\nParticle 1 initial: pos={system.positions[0]}, vel={system.velocities[0]}")
    print(f"Particle 2 initial: pos={system.positions[1]}, vel={system.velocities[1]}")

    # Record trajectory
    times = []
    pos_x1 = []
    pos_x2 = []
    vel_x1 = []
    vel_x2 = []

    t_end = 2.0
    n_steps = int(t_end / system.dt)

    for _ in range(n_steps):
        times.append(system.time)
        pos_x1.append(system.positions[0, 0])
        pos_x2.append(system.positions[1, 0])
        vel_x1.append(system.velocities[0, 0])
        vel_x2.append(system.velocities[1, 0])
        system.step()

    print(f"\nParticle 1 final: pos={system.positions[0]}, vel={system.velocities[0]}")
    print(f"Particle 2 final: pos={system.positions[1]}, vel={system.velocities[1]}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(times, pos_x1, "b-", linewidth=2, label="Particle 1")
    ax1.plot(times, pos_x2, "r-", linewidth=2, label="Particle 2")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("X Position [m]")
    ax1.set_title("Particle Positions")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(times, vel_x1, "b-", linewidth=2, label="Particle 1")
    ax2.plot(times, vel_x2, "r-", linewidth=2, label="Particle 2")
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("X Velocity [m/s]")
    ax2.set_title("Particle Velocities")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("dem_particle_collision.png", dpi=150, bbox_inches="tight")
    print("Plot saved as: dem_particle_collision.png")

    plt.show()


def verify_energy_conservation():
    """Verify energy dissipation in the system."""
    print("\n" + "=" * 60)
    print("DEM Verification - Energy Analysis")
    print("=" * 60)

    system = ParticleSystem(
        n_particles=30,
        domain_size=(5.0, 8.0),
        particle_radius=0.15,
        particle_density=2500,
        gravity=9.81,
        k_n=1e5,
        damping=0.3,
        dt=1e-4,
    )

    print(f"\n{system.n_particles} particles in {system.width}x{system.height}m domain")

    # Record energy over time
    times = []
    kinetic_energies = []

    t_end = 3.0
    n_steps = int(t_end / system.dt)
    record_interval = 100

    for i in range(n_steps):
        if i % record_interval == 0:
            times.append(system.time)
            kinetic_energies.append(system.compute_kinetic_energy())
        system.step()

    print(f"\nInitial KE: {kinetic_energies[0]:.2e}J")
    print(f"Final KE: {kinetic_energies[-1]:.2e}J")
    print(f"Energy dissipation: {(1 - kinetic_energies[-1]/kinetic_energies[0])*100:.1f}%")

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, kinetic_energies, "b-", linewidth=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Kinetic Energy [J]")
    ax.set_title("System Kinetic Energy vs Time")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("dem_energy_dissipation.png", dpi=150, bbox_inches="tight")
    print("\nPlot saved as: dem_energy_dissipation.png")

    plt.show()


def main():
    """Run all verification tests."""
    print("\n🌊 DEM Solver - Verification Suite\n")

    # Test 1: Single particle fall
    verify_single_particle_fall()

    # Test 2: Particle collision
    verify_particle_collision()

    # Test 3: Energy dissipation
    verify_energy_conservation()

    print("\n" + "=" * 60)
    print("✅ Verification Complete!")
    print("=" * 60)
    print("\nNotes:")
    print("- Single particle should accelerate under gravity until hitting bottom")
    print("- Two particles should bounce off each other due to contact forces")
    print("- Energy should dissipate over time due to damping")
    print("- This is a simplified DEM for educational purposes")


if __name__ == "__main__":
    main()
