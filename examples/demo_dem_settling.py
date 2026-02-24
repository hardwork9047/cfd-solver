"""
Demo: Particle Settling Simulation using DEM

This demo shows particles falling and settling under gravity.
"""

import matplotlib.pyplot as plt
import os
import shutil

from dem import ParticleSystem


def demo_particle_settling():
    """Demonstrate particle settling under gravity."""
    print("🌊 DEM Solver - Particle Settling Demo")
    print("=" * 60)
    print("Simulating particle settling under gravity...")
    print()

    # Create particle system with fewer particles for stability
    system = ParticleSystem(
        n_particles=4,
        domain_size=(2.0, 3.0),
        particle_radius=0.1,
        particle_density=2500,
        gravity=9.81,
        k_n=1e4,
        damping=0.8,
        dt=1e-4,
    )

    print(f"Number of particles: {system.n_particles}")
    print(f"Domain size: {system.width}m × {system.height}m")
    print(f"Particle radius: {system.radius}m")
    print(f"Particle mass: {system.mass:.3e}kg")
    print()

    # Set initial positions manually for better visualization
    system.positions[0] = [0.5, 2.5]
    system.positions[1] = [1.5, 2.7]
    system.positions[2] = [0.8, 2.0]
    system.positions[3] = [1.2, 2.3]
    
    # Create output directories
    results_dir = "results"
    frames_dir = os.path.join(results_dir, "frames")
    images_dir = os.path.join(results_dir, "images")
    video_dir = os.path.join(results_dir, "video")
    
    for directory in [results_dir, frames_dir, images_dir, video_dir]:
        os.makedirs(directory, exist_ok=True)
    
    # Clean up old frames
    if os.path.exists(frames_dir):
        for file in os.listdir(frames_dir):
            if file.startswith("dem_frame_"):
                os.remove(os.path.join(frames_dir, file))
    
    # Plot initial state
    print("Initial state:")
    fig_init = system.plot(save_path=os.path.join(frames_dir, "dem_frame_0000.png"))
    plt.close(fig_init)

    # Run simulation with frame capture
    print("\nRunning simulation and generating frames...")
    t_end = 2.0
    frame_interval = 200  # Save frame every 200 iterations
    frame_count = 0
    
    n_steps = int(t_end / system.dt)
    
    for step in range(n_steps):
        system.step()
        
        if step % frame_interval == 0:
            frame_count += 1
            fig = system.plot(save_path=os.path.join(frames_dir, f"dem_frame_{frame_count:04d}.png"))
            plt.close(fig)
            
        if step % 2000 == 0:
            ke = system.compute_kinetic_energy()
            print(f"  Step {step}/{n_steps}, Time {system.time:.3f}s, KE={ke:.2e}J")

    # Final frame
    frame_count += 1
    fig_final = system.plot(save_path=os.path.join(frames_dir, f"dem_frame_{frame_count:04d}.png"))
    plt.close(fig_final)

    print(f"\n✅ Generated {frame_count + 1} frames in '{frames_dir}/' directory")
    
    video_path = os.path.join(video_dir, "dem_settling.mp4")
    print("\nTo create video, run:")
    print(f"  ffmpeg -framerate 10 -i {frames_dir}/dem_frame_%04d.png -c:v libx264 -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" {video_path} -y")
    
    # Plot particle distribution
    bin_centers, hist = system.get_particle_distribution(n_bins=10)

    fig_dist, ax = plt.subplots(figsize=(8, 6))
    ax.barh(bin_centers, hist, height=bin_centers[1] - bin_centers[0], alpha=0.7, color='steelblue')
    ax.set_xlabel("Number of particles")
    ax.set_ylabel("Height [m]")
    ax.set_title(f"Vertical Distribution at t={system.time:.2f}s")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    dist_path = os.path.join(images_dir, "dem_distribution.png")
    plt.savefig(dist_path, dpi=150, bbox_inches="tight")
    print(f"\nDistribution plot saved as: {dist_path}")

    plt.show()

    print("\n✅ Demo complete!")
    print("\nOutput locations:")
    print(f"  - Frames: {frames_dir}/")
    print(f"  - Images: {images_dir}/")
    print(f"  - Video: {video_dir}/")
    print("\nWhat to observe in the video:")
    print("- 4 particles initially positioned at different heights")
    print("- Particles fall under gravity")
    print("- Particles collide with each other and walls")
    print("- Eventually settle at the bottom")
    print("- Energy dissipates due to damping")


if __name__ == "__main__":
    demo_particle_settling()
