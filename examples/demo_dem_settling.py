"""
Demo: Particle Settling Simulation using DEM

This demo shows particles falling and settling under gravity.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
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
        domain_size=(0.6, 0.8),
        particle_radius=0.1,
        particle_density=2500,
        gravity=9.81,
        k_n=1e6,  # Very high stiffness to minimize penetration
        damping=0.8,
        dt=1e-5,
    )

    print(f"Number of particles: {system.n_particles}")
    print(f"Domain size: {system.width}m × {system.height}m")
    print(f"Particle radius: {system.radius}m")
    print(f"Particle mass: {system.mass:.3e}kg")
    print()

    # Set initial positions without overlap - manual placement for 4 particles
    # Place particles in an unbalanced configuration to see rolling behavior
    # Distance between particles must be > 0.2m (diameter)
    system.positions[0] = [0.18, 0.12]  # Bottom - left side
    system.positions[1] = [0.42, 0.12]  # Bottom - right side (asymmetric)
    system.positions[2] = [0.25, 0.45]  # Top left - offset from bottom particles
    system.positions[3] = [0.38, 0.60]  # Top right - higher and offset
    
    # Verify no initial overlaps
    print("\nInitial particle positions:")
    for i in range(system.n_particles):
        print(f"  Particle {i}: {system.positions[i]}")
    
    print("\nInitial distances:")
    for i in range(system.n_particles):
        for j in range(i+1, system.n_particles):
            dist = np.linalg.norm(system.positions[j] - system.positions[i])
            overlap = 2*system.radius - dist
            status = "OVERLAP!" if overlap > 0 else "OK"
            print(f"  Particles {i}-{j}: distance={dist:.4f}m, overlap={overlap:.4f}m [{status}]")
    
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

    # Open CSV file for logging particle data
    csv_path = os.path.join(results_dir, "particle_data.csv")
    csv_file = open(csv_path, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    
    # Write CSV header
    header = ['step', 'time']
    for i in range(system.n_particles):
        header.extend([f'p{i}_x', f'p{i}_y', f'p{i}_vx', f'p{i}_vy', f'p{i}_fx', f'p{i}_fy'])
    csv_writer.writerow(header)
    
    # Write initial state
    row = [0, 0.0]
    for i in range(system.n_particles):
        row.extend([
            system.positions[i, 0], system.positions[i, 1],
            system.velocities[i, 0], system.velocities[i, 1],
            0.0, -system.mass * system.gravity  # Initial forces (just gravity)
        ])
    csv_writer.writerow(row)

    # Run simulation with frame capture
    print("\nRunning simulation and generating frames...")
    t_end = 0.5
    frame_interval = 200  # Save frame every 200 iterations
    frame_count = 0
    
    n_steps = int(t_end / system.dt)
    
    for step in range(n_steps):
        system.step()
        
        # Write particle data to CSV
        row = [step + 1, system.time]
        for i in range(system.n_particles):
            row.extend([
                system.positions[i, 0], system.positions[i, 1],
                system.velocities[i, 0], system.velocities[i, 1],
                system.forces[i, 0], system.forces[i, 1]
            ])
        csv_writer.writerow(row)
        
        if step % frame_interval == 0:
            frame_count += 1
            fig = system.plot(save_path=os.path.join(frames_dir, f"dem_frame_{frame_count:04d}.png"))
            plt.close(fig)
            
        if step % 5000 == 0:
            ke = system.compute_kinetic_energy()
            print(f"  Step {step}/{n_steps}, Time {system.time:.3f}s, KE={ke:.2e}J")
    
    # Close CSV file
    csv_file.close()
    print(f"\nParticle data saved to: {csv_path}")

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
    print("- 20 particles initially distributed without overlap in upper half")
    print("- Particles fall under gravity")
    print("- Particles collide with each other and walls")
    print("- Color indicates force magnitude (blue=low, red=high)")
    print("- Eventually settle at the bottom")
    print("- Energy dissipates due to damping")


if __name__ == "__main__":
    demo_particle_settling()
