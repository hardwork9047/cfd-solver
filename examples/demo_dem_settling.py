"""
Demo: Particle Settling Simulation using DEM

This demo shows particles falling and settling under gravity.
"""

import matplotlib.pyplot as plt

from dem import ParticleSystem


def demo_particle_settling():
    """Demonstrate particle settling under gravity."""
    print("🌊 DEM Solver - Particle Settling Demo")
    print("=" * 60)
    print("Simulating particle settling under gravity...")
    print()

    # Create particle system
    system = ParticleSystem(
        n_particles=50,
        domain_size=(5.0, 10.0),
        particle_radius=0.15,
        particle_density=2500,
        gravity=9.81,
        k_n=1e5,
        damping=0.5,
        dt=1e-4,
    )

    print(f"Number of particles: {system.n_particles}")
    print(f"Domain size: {system.width}m × {system.height}m")
    print(f"Particle radius: {system.radius}m")
    print(f"Particle mass: {system.mass:.3e}kg")
    print()

    # Plot initial state
    print("Initial state:")
    fig_init = system.plot(save_path="dem_initial.png")
    plt.close(fig_init)

    # Run simulation
    print("\nRunning simulation...")
    system.run(t_end=5.0, output_interval=2000)

    # Plot final state
    print("\nFinal state:")
    fig_final = system.plot(save_path="dem_final.png")

    # Plot particle distribution
    bin_centers, hist = system.get_particle_distribution(n_bins=20)

    fig_dist, ax = plt.subplots(figsize=(8, 6))
    ax.barh(bin_centers, hist, height=bin_centers[1] - bin_centers[0], alpha=0.7)
    ax.set_xlabel("Number of particles")
    ax.set_ylabel("Height [m]")
    ax.set_title("Vertical Distribution of Particles")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("dem_distribution.png", dpi=150, bbox_inches="tight")
    print("Distribution plot saved as: dem_distribution.png")

    # Show plots
    plt.show()

    print("\n✅ Demo complete!")
    print("\nWhat to observe:")
    print("- Particles initially distributed randomly fall under gravity")
    print("- Particles collide with each other and walls")
    print("- Eventually, particles settle at the bottom forming a packed bed")
    print("- Energy dissipates due to damping")


if __name__ == "__main__":
    demo_particle_settling()
