"""
Demo: Large-scale Particle Settling Simulation with 100 particles

This demo shows 100 particles falling and settling under gravity in a larger domain
with 5x longer simulation time.
粒子径を乱数で±15%変化させています。
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

from dem import ParticleSystem


def demo_large_particle_settling():
    """Demonstrate particle settling with 100 particles in a larger domain."""
    print("🌊 DEM Solver - Large-scale Particle Settling Demo (100 particles)")
    print("=" * 60)
    print("Simulating 100 particles falling and settling under gravity...")
    print("粒子径：基本値 ±15% の範囲で乱数変化")
    print()

    # Create particle system with 100 particles in a larger domain
    system = ParticleSystem(
        n_particles=100,
        domain_size=(2.5, 4.0),
        particle_radius=0.08,
        particle_density=2500,
        gravity=9.81,
        k_n=5e6,  # Very high stiffness to minimize penetration
        damping=0.8,
        dt=1e-5,
    )
    
    # 粒子径を乱数で変化させる（±15%）
    base_radius = 0.08
    diameter_variation = 0.15
    min_radius = base_radius * (1 - diameter_variation)
    max_radius = base_radius * (1 + diameter_variation)
    
    np.random.seed(42)  # For reproducibility
    radii = np.random.uniform(min_radius, max_radius, system.n_particles)
    system.radii = radii
    
    # 粒子径に応じて質量を更新
    system.masses = system.density * np.pi * radii**2
    
    print(f"粒子径統計:")
    print(f"  基本粒子径: {base_radius:.4f} m")
    print(f"  粒子径範囲: {min_radius:.4f} - {max_radius:.4f} m")
    print(f"  平均粒子径: {np.mean(radii):.4f} m")
    print(f"  標準偏差: {np.std(radii):.4f} m")
    print(f"  最小粒子径: {np.min(radii):.4f} m")
    print(f"  最大粒子径: {np.max(radii):.4f} m")
    print()
    
    # 設定確認
    print("設定確認:")
    print(f"  system.radii の長さ: {len(system.radii)}")
    print(f"  system.radii の最初の5個: {system.radii[:5]}")
    print(f"  system.masses の最初の5個: {system.masses[:5]}")
    print()

    print(f"Number of particles: {system.n_particles}")
    print(f"Domain size: {system.width}m × {system.height}m")
    print()

    # Set initial positions without overlap - random placement in upper half
    print("初期位置を設定中（重なり判定あり）...")
    placed = 0
    max_total_attempts = 100000
    global_attempt = 0
    
    for i in range(system.n_particles):
        placed_i = False
        max_attempts = 10000
        
        for attempt in range(max_attempts):
            global_attempt += 1
            if global_attempt > max_total_attempts:
                print(f"警告: 最大試行回数に達しました。{placed}/{system.n_particles} 個の粒子を配置しました。")
                break
            
            # Generate random position in upper half with margin for radius
            x = np.random.uniform(radii[i], system.width - radii[i])
            y = np.random.uniform(system.height * 0.5, system.height - radii[i])
            
            # Check distance to all existing particles (using variable radii)
            valid = True
            for j in range(i):
                dx = x - system.positions[j, 0]
                dy = y - system.positions[j, 1]
                dist = np.sqrt(dx**2 + dy**2)
                min_distance = radii[i] + radii[j] + 1e-6  # 重なりなし + 余裕
                if dist < min_distance:
                    valid = False
                    break
            
            if valid:
                system.positions[i] = [x, y]
                placed_i = True
                placed += 1
                break
        
        if not placed_i:
            print(f"警告: 粒子 {i} の配置に失敗しました（{max_attempts}回の試行後）")
        
        if (i + 1) % 20 == 0:
            print(f"  {i + 1}/{system.n_particles} 個の粒子を配置完了")
    
    print(f"✓ {placed}/{system.n_particles} 個の粒子を配置完了\n")
    
    # Verify no initial overlaps
    print("重なり判定を実行中...")
    overlap_count = 0
    max_overlap = 0.0
    
    for i in range(system.n_particles):
        for j in range(i+1, system.n_particles):
            dx = system.positions[j, 0] - system.positions[i, 0]
            dy = system.positions[j, 1] - system.positions[i, 1]
            dist = np.sqrt(dx**2 + dy**2)
            min_allowed_dist = radii[i] + radii[j]
            overlap = min_allowed_dist - dist
            
            if overlap > 1e-6:  # 誤差範囲以上の重なり
                overlap_count += 1
                max_overlap = max(max_overlap, overlap)
                if overlap_count <= 10:  # 最初の10個だけ表示
                    print(f"  警告: 粒子 {i}-{j}: 距離={dist:.6f}m, 重なり={overlap:.6f}m")
    
    if overlap_count == 0:
        print("✓ すべての粒子が重なりなく初期化されました")
    else:
        print(f"⚠️  {overlap_count} 個の重なるペアが検出されました")
        print(f"   最大重なり量: {max_overlap:.6f}m")
    print()
    
    # Create output directories
    results_dir = "results"
    frames_dir = os.path.join(results_dir, "frames_100p")
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
    csv_path = os.path.join(results_dir, "particle_data_100p.csv")
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

    # Run simulation with frame capture (5x longer simulation time)
    print("\nRunning simulation and generating frames...")
    t_end = 5.0  # 5x longer than 50-particle demo
    frame_interval = 2000  # Save frame every 2000 iterations
    csv_interval = 500  # Write CSV every 500 iterations
    frame_count = 0
    
    n_steps = int(t_end / system.dt)
    
    print(f"Total timesteps: {n_steps}")
    print()
    
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
            
        if step % 20000 == 0:
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
    
    video_path = os.path.join(video_dir, "dem_settling_100p.mp4")
    print("\nTo create video, run:")
    print(f"  ffmpeg -framerate 10 -i {frames_dir}/dem_frame_%04d.png -c:v libx264 -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" {video_path} -y")
    
    # Plot particle distribution
    bin_centers, hist = system.get_particle_distribution(n_bins=20)

    fig_dist, ax = plt.subplots(figsize=(8, 6))
    ax.barh(bin_centers, hist, height=bin_centers[1] - bin_centers[0], alpha=0.7, color='steelblue')
    ax.set_xlabel("Number of particles")
    ax.set_ylabel("Height [m]")
    ax.set_title(f"Vertical Distribution at t={system.time:.2f}s (100 particles)")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    dist_path = os.path.join(images_dir, "dem_distribution_100p.png")
    plt.savefig(dist_path, dpi=150, bbox_inches="tight")
    print(f"\nDistribution plot saved as: {dist_path}")

    plt.show()

    print("\n✅ Demo complete!")
    print("\nOutput locations:")
    print(f"  - Frames: {frames_dir}/")
    print(f"  - Images: {images_dir}/")
    print(f"  - Video: {video_dir}/")
    print("\nSimulation parameters:")
    print(f"  - Particles: {system.n_particles}")
    print(f"  - Total time: {t_end}s (5x longer)")
    print(f"  - Total steps: {n_steps:,}")
    print("\nWhat to observe in the video:")
    print("- 100 particles initially distributed without overlap in upper half")
    print("- Particles fall under gravity")
    print("- Particles collide with each other and walls")
    print("- Color indicates force magnitude (blue=low, red=high)")
    print("- Eventually settle at the bottom forming a packed bed")
    print("- Energy dissipates due to damping")


if __name__ == "__main__":
    demo_large_particle_settling()
