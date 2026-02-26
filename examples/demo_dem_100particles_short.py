"""
Demo: Large-scale Particle Settling Simulation with 100 particles (SHORT VERSION)

This is a shortened version for testing purposes.
粒子径を乱数で±15%変化させています。
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

from dem import ParticleSystem


def demo_large_particle_settling_short():
    """Demonstrate particle settling with 100 particles - SHORT VERSION."""
    print("🌊 DEM Solver - Large-scale Particle Settling Demo (100 particles) - SHORT VERSION")
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
        k_n=5e6,
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

    print(f"Number of particles: {system.n_particles}")
    print(f"Domain size: {system.width}m × {system.height}m")
    print()

    # Set initial positions without overlap
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
            
            x = np.random.uniform(radii[i], system.width - radii[i])
            y = np.random.uniform(system.height * 0.5, system.height - radii[i])
            
            valid = True
            for j in range(i):
                dx = x - system.positions[j, 0]
                dy = y - system.positions[j, 1]
                dist = np.sqrt(dx**2 + dy**2)
                min_distance = radii[i] + radii[j] + 1e-6
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
            
            if overlap > 1e-6:
                overlap_count += 1
                max_overlap = max(max_overlap, overlap)
                if overlap_count <= 10:
                    print(f"  警告: 粒子 {i}-{j}: 距離={dist:.6f}m, 重なり={overlap:.6f}m")
    
    if overlap_count == 0:
        print("✓ すべての粒子が重なりなく初期化されました")
    else:
        print(f"⚠️  {overlap_count} 個の重なるペアが検出されました")
        print(f"   最大重なり量: {max_overlap:.6f}m")
    print()
    
    # Create output directories
    results_dir = "results"
    frames_dir = os.path.join(results_dir, "frames_100p_short")
    images_dir = os.path.join(results_dir, "images")
    
    for directory in [results_dir, frames_dir, images_dir]:
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

    # Run SHORT simulation
    print("\nRunning SHORT simulation (10 seconds, 1 million steps)...")
    t_end = 0.5  # SHORT: 0.5秒のみ（フル版は5秒）
    frame_interval = 25000  # Save frame every 25000 iterations
    frame_count = 0
    
    n_steps = int(t_end / system.dt)
    print(f"Total timesteps: {n_steps:,}")
    print()
    
    for step in range(n_steps):
        system.step()
        
        if step % frame_interval == 0:
            frame_count += 1
            print(f"  Step {step:,}/{n_steps:,}, Time {system.time:.3f}s")
            fig = system.plot(save_path=os.path.join(frames_dir, f"dem_frame_{frame_count:04d}.png"))
            plt.close(fig)
        
        if step % 100000 == 0 and step > 0:
            ke = system.compute_kinetic_energy()
            print(f"    KE={ke:.2e}J")
    
    # Final frame
    frame_count += 1
    fig_final = system.plot(save_path=os.path.join(frames_dir, f"dem_frame_{frame_count:04d}.png"))
    plt.close(fig_final)

    print(f"\n✅ Generated {frame_count + 1} frames in '{frames_dir}/' directory")
    
    print("\n✅ Demo complete!")
    print("\nOutput locations:")
    print(f"  - Frames: {frames_dir}/")
    print(f"  - Images: {images_dir}/")
    print("\nSimulation parameters:")
    print(f"  - Particles: {system.n_particles}")
    print(f"  - Particle radii variation: ±{diameter_variation*100:.1f}%")
    print(f"  - Total time: {t_end}s (SHORT VERSION)")
    print(f"  - Total steps: {n_steps:,}")
    print(f"  - Grid optimization: 空間パーティショニング使用")


if __name__ == "__main__":
    demo_large_particle_settling_short()
