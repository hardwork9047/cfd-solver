"""
Demo: Large-scale Particle Settling Simulation with 50 particles

This demo shows 50 particles falling and settling under gravity in a larger domain.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

from dem import ParticleSystem


def demo_large_particle_settling():
    """Demonstrate particle settling with 50 particles in a larger domain."""
    print("🌊 DEM Solver - Large-scale Particle Settling Demo (50 particles)")
    print("=" * 60)
    print("Simulating 50 particles falling and settling under gravity...")
    print()

    # Create particle system with 50 particles in a larger domain
    system = ParticleSystem(
        n_particles=50,
        domain_size=(2.0, 3.0),
        particle_radius=0.08,
        particle_density=2500,
        gravity=9.81,
        k_n=5e6,  # Very high stiffness to minimize penetration
        damping=0.8,
        dt=1e-5,
    )

    print(f"Number of particles: {system.n_particles}")
    print(f"Domain size: {system.width}m × {system.height}m")
    print(f"Particle radius: {system.radius}m")
    print(f"Particle mass: {system.mass:.3e}kg")
    print()

    # Set initial positions without overlap - random placement in upper half
    np.random.seed(42)  # For reproducibility
    min_distance = 2.2 * system.radius  # Minimum distance between particle centers
    
    for i in range(system.n_particles):
        max_attempts = 1000
        for attempt in range(max_attempts):
            # Generate random position in upper half
            x = np.random.uniform(system.radius, system.width - system.radius)
            y = np.random.uniform(system.height * 0.5, system.height - system.radius)
            
            # Check distance to all existing particles
            valid = True
            for j in range(i):
                dx = x - system.positions[j, 0]
                dy = y - system.positions[j, 1]
                dist = np.sqrt(dx**2 + dy**2)
                if dist < min_distance:
                    valid = False
                    break
            
            if valid:
                system.positions[i] = [x, y]
                break
        else:
            # If no valid position found, place it anyway (fallback)
            system.positions[i] = [x, y]
            print(f"Warning: Could not find non-overlapping position for particle {i}")
    
    # Verify no initial overlaps
    print("\nChecking initial overlaps...")
    overlap_count = 0
    for i in range(system.n_particles):
        for j in range(i+1, system.n_particles):
            dist = np.linalg.norm(system.positions[j] - system.positions[i])
            overlap = 2*system.radius - dist
            if overlap > 0:
                overlap_count += 1
                print(f"  WARNING: Particles {i}-{j}: distance={dist:.4f}m, overlap={overlap:.4f}m")
    
    if overlap_count == 0:
        print("✓ All particles initialized without overlap")
    else:
        print(f"⚠️  Found {overlap_count} overlapping pairs")
    print()
    
    # Create output directories
    results_dir = "results"
    frames_dir = os.path.join(results_dir, "frames_50p")
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
    print("Saving initial state...")
    fig_init = system.plot(save_path=os.path.join(frames_dir, "dem_frame_0000.png"))
    plt.close(fig_init)

    # Open CSV file for logging particle data (limited to save disk space)
    csv_path = os.path.join(results_dir, "particle_data_50p.csv")
    csv_file = open(csv_path, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    
    # Write CSV header (only positions and velocities to save space)
    header = ['step', 'time', 'kinetic_energy']
    for i in range(min(10, system.n_particles)):  # Only log first 10 particles
        header.extend([f'p{i}_x', f'p{i}_y', f'p{i}_vx', f'p{i}_vy'])
    csv_writer.writerow(header)
    
    # Write initial state
    row = [0, 0.0, 0.0]
    for i in range(min(10, system.n_particles)):
        row.extend([
            system.positions[i, 0], system.positions[i, 1],
            system.velocities[i, 0], system.velocities[i, 1]
        ])
    csv_writer.writerow(row)

    # Run simulation with frame capture
    print("\nRunning simulation and generating frames...")
    t_end = 1.0
    frame_interval = 400  # Save frame every 400 iterations
    csv_interval = 100  # Write CSV every 100 iterations
    frame_count = 0
    
    n_steps = int(t_end / system.dt)
    
    for step in range(n_steps):
        system.step()
        
        # Write particle data to CSV (periodically)
        if step % csv_interval == 0:
            ke = system.compute_kinetic_energy()
            row = [step + 1, system.time, ke]
            for i in range(min(10, system.n_particles)):
                row.extend([
                    system.positions[i, 0], system.positions[i, 1],
                    system.velocities[i, 0], system.velocities[i, 1]
                ])
            csv_writer.writerow(row)
        
        if step % frame_interval == 0:
            frame_count += 1
            fig = system.plot(save_path=os.path.join(frames_dir, f"dem_frame_{frame_count:04d}.png"))
            plt.close(fig)
            
        if step % 10000 == 0:
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
    
    video_path = os.path.join(video_dir, "dem_settling_50p.mp4")
    print("\nTo create video, run:")
    print(f"  ffmpeg -framerate 10 -i {frames_dir}/dem_frame_%04d.png -c:v libx264 -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" {video_path} -y")
    
    # Plot particle distribution
    bin_centers, hist = system.get_particle_distribution(n_bins=15)

    fig_dist, ax = plt.subplots(figsize=(8, 6))
    ax.barh(bin_centers, hist, height=bin_centers[1] - bin_centers[0], alpha=0.7, color='steelblue')
    ax.set_xlabel("Number of particles")
    ax.set_ylabel("Height [m]")
    ax.set_title(f"Vertical Distribution at t={system.time:.2f}s (50 particles)")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    dist_path = os.path.join(images_dir, "dem_distribution_50p.png")
    plt.savefig(dist_path, dpi=150, bbox_inches="tight")
    print(f"\nDistribution plot saved as: {dist_path}")

    plt.show()

    print("\n✅ Demo complete!")
    print("\nOutput locations:")
    print(f"  - Frames: {frames_dir}/")
    print(f"  - Images: {images_dir}/")
    print(f"  - Video: {video_dir}/")
    print("\nWhat to observe in the video:")
    print("- 50 particles initially distributed without overlap in upper half")
    print("- Particles fall under gravity")
    print("- Particles collide with each other and walls")
    print("- Color indicates force magnitude (blue=low, red=high)")
    print("- Eventually settle at the bottom forming a packed bed")
    print("- Energy dissipates due to damping")


if __name__ == "__main__":
    demo_large_particle_settling()
