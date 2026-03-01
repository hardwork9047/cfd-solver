#!/usr/bin/env python
"""
粒子シミュレーション：Stokes抗力 + 粒子径乱数

物理モデル:
  - 重力 + Stokes抗力 (F_drag = -3πμDv)
  - 陰解法による安定な時間積分
  - 終端速度: v_term = ρ_p * D² * g / (18μ_f)  → 粒子径に強依存

パラメータ：
- 粒子数: 500個
- 粒子径: 基本値 ±50% の範囲で乱数変化
- タイムステップ数: 1500ステップ (1.5s)
"""

import logging
import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


class ParticleSimulation:
    """Stokes抗力を持つ粒子シミュレーション"""

    def __init__(
        self,
        num_particles=500,
        base_diameter=1.0e-6,       # 基本粒子径 [m]
        diameter_variation=0.50,    # 粒子径変動 ±50%
        domain_size=1.0e-4,         # シミュレーション領域 [m]
        timesteps=1500,
        dt=1.0e-3,                  # 時間ステップ [s]（陰解法で大きくても安定）
        mu_fluid=1.8e-5,            # 空気の粘性係数 [Pa·s]
        rho_fluid=1.2,              # 空気密度 [kg/m³]
        rho_particle=1000,          # 粒子密度 [kg/m³]（水）
    ):
        self.num_particles = num_particles
        self.base_diameter = base_diameter
        self.diameter_variation = diameter_variation
        self.domain_size = domain_size
        self.timesteps = timesteps
        self.dt = dt
        self.mu_fluid = mu_fluid
        self.rho_fluid = rho_fluid
        self.rho_particle = rho_particle

        # 粒子座標・速度
        self.positions = np.random.uniform(0, domain_size, size=(num_particles, 2))
        self.velocities = np.random.normal(0, 0.01, size=(num_particles, 2))

        # 粒子径（±50%）
        min_d = base_diameter * (1 - diameter_variation)
        max_d = base_diameter * (1 + diameter_variation)
        self.diameters = np.random.uniform(min_d, max_d, size=num_particles)

        # 質量
        volume = (4 / 3) * np.pi * (self.diameters / 2) ** 3
        self.masses = rho_particle * volume

        # Stokes緩和時間: τ = ρ_p * D² / (18 * μ_f)
        # 大きい粒子ほど長い → 終端速度に到達するのが遅い
        self.taus = rho_particle * self.diameters**2 / (18 * mu_fluid)

        # 終端速度: v_term = τ * g  (D² に比例 → 大きな粒子ほど速く落下)
        g_vec = np.array([0.0, -9.8])
        self.v_terminal = np.outer(self.taus, g_vec)  # shape (N, 2)

        self.time = 0.0
        self.step_count = 0

        logger.info("粒子シミュレーション初期化:")
        logger.info(f"  粒子数: {num_particles}")
        logger.info(f"  粒子径範囲: {min_d*1e6:.2f} – {max_d*1e6:.2f} μm")
        logger.info(f"  τ範囲: {np.min(self.taus)*1e6:.1f} – {np.max(self.taus)*1e6:.1f} μs")
        vt_min = np.min(np.abs(self.v_terminal[:, 1])) * 1e6
        vt_max = np.max(np.abs(self.v_terminal[:, 1])) * 1e6
        logger.info(f"  終端速度範囲: {vt_min:.1f} – {vt_max:.1f} μm/s")

    # ------------------------------------------------------------------
    def step(self):
        """1ステップ: Stokes抗力の陰解法 (線形抵抗の厳密解) + 壁反射"""
        dt = self.dt
        # 厳密解: v(t+dt) = v_term + (v - v_term) * exp(-dt/τ)
        decay = np.exp(-dt / self.taus)[:, None]     # shape (N, 1)
        v_excess = self.velocities - self.v_terminal  # v - v_term

        # 位置の厳密積分: Δx = v_term*dt + v_excess*τ*(1 - exp(-dt/τ))
        self.positions += (
            self.v_terminal * dt
            + v_excess * self.taus[:, None] * (1.0 - decay)
        )
        self.velocities = self.v_terminal + v_excess * decay

        # 壁反射（反発係数 0.8）
        for dim in range(2):
            lo = self.positions[:, dim] < 0
            hi = self.positions[:, dim] > self.domain_size
            self.positions[lo, dim] = 0.0
            self.positions[hi, dim] = self.domain_size
            self.velocities[lo, dim] *= -0.8
            self.velocities[hi, dim] *= -0.8

        self.time += dt
        self.step_count += 1

    # ------------------------------------------------------------------
    def get_drag_forces(self):
        """各粒子の Stokes 抗力の大きさ [N] を返す

        F_drag = 3π * μ_f * D * |v|
        終端速度では F_drag = m*g ∝ D³ → 粒子径が大きいほど抗力が大きい
        """
        v_mag = np.linalg.norm(self.velocities, axis=1)
        return 3 * np.pi * self.mu_fluid * self.diameters * v_mag

    # ------------------------------------------------------------------
    def run(self):
        """シミュレーション全ステップを実行"""
        logger.info("シミュレーション実行中...")
        for s in range(self.timesteps):
            self.step()
            if (s + 1) % 300 == 0:
                logger.info(
                    f"  step {s+1}/{self.timesteps}  t={self.time:.4e}s"
                )
        logger.info("完了")

    # ------------------------------------------------------------------
    def visualize(self):
        """粒子を Stokes 抗力で色付けして可視化"""
        drag_fN = self.get_drag_forces() * 1e15  # fN (フェムトニュートン)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # --- 左: 抗力で色付け ---
        sc1 = axes[0].scatter(
            self.positions[:, 0] * 1e6,
            self.positions[:, 1] * 1e6,
            s=30 * (self.diameters / self.base_diameter) ** 2,
            c=drag_fN,
            cmap="plasma",
            vmin=0,
            vmax=np.percentile(drag_fN, 98),
            alpha=0.75,
            edgecolors="none",
        )
        axes[0].set_xlim(0, self.domain_size * 1e6)
        axes[0].set_ylim(0, self.domain_size * 1e6)
        axes[0].set_aspect("equal")
        axes[0].set_xlabel("X [μm]"); axes[0].set_ylabel("Y [μm]")
        axes[0].set_title(
            f"Colored by Stokes Drag  t={self.time:.3f}s  N={self.num_particles}"
        )
        cb1 = plt.colorbar(sc1, ax=axes[0])
        cb1.set_label("Stokes drag [fN]")

        # --- right: colored by diameter (reference) ---
        sc2 = axes[1].scatter(
            self.positions[:, 0] * 1e6,
            self.positions[:, 1] * 1e6,
            s=30 * (self.diameters / self.base_diameter) ** 2,
            c=self.diameters * 1e6,
            cmap="viridis",
            alpha=0.75,
            edgecolors="none",
        )
        axes[1].set_xlim(0, self.domain_size * 1e6)
        axes[1].set_ylim(0, self.domain_size * 1e6)
        axes[1].set_aspect("equal")
        axes[1].set_xlabel("X [μm]"); axes[1].set_ylabel("Y [μm]")
        axes[1].set_title("Colored by Diameter (reference)")
        cb2 = plt.colorbar(sc2, ax=axes[1])
        cb2.set_label("D [μm]")

        plt.suptitle(
            f"Particle Simulation  Re_p<<1 (Stokes regime)  "
            f"D={self.base_diameter*1e6:.0f}um +/-{self.diameter_variation*100:.0f}%  "
            f"step={self.step_count}",
            fontsize=11,
        )
        plt.tight_layout()
        return fig

    # ------------------------------------------------------------------
    def get_statistics(self):
        """粒子の統計情報"""
        v_mag = np.linalg.norm(self.velocities, axis=1)
        drag = self.get_drag_forces()
        return {
            "mean_diameter": np.mean(self.diameters) * 1e6,
            "std_diameter":  np.std(self.diameters)  * 1e6,
            "min_diameter":  np.min(self.diameters)  * 1e6,
            "max_diameter":  np.max(self.diameters)  * 1e6,
            "mean_velocity": np.mean(v_mag),
            "max_velocity":  np.max(v_mag),
            "mean_drag_fN":  np.mean(drag) * 1e15,
            "max_drag_fN":   np.max(drag)  * 1e15,
            "kinetic_energy": 0.5 * np.sum(self.masses * v_mag**2),
        }


# ======================================================================
def main():
    print("\n" + "=" * 70)
    print("粒子シミュレーション：Stokes抗力 + 粒子径乱数")
    print("=" * 70)

    sim = ParticleSimulation(
        num_particles=500,
        base_diameter=1.0e-6,
        diameter_variation=0.50,
        domain_size=1.0e-4,
        timesteps=1500,
        dt=1.0e-3,
    )

    # 初期統計
    stats0 = sim.get_statistics()
    print(f"\n初期状態:")
    print(f"  粒子径: {stats0['min_diameter']:.2f} – {stats0['max_diameter']:.2f} μm  "
          f"(mean={stats0['mean_diameter']:.2f} μm)")
    print(f"  抗力:   {stats0['mean_drag_fN']:.4f} fN (mean)  "
          f"{stats0['max_drag_fN']:.4f} fN (max)")
    print(f"  ※ 抗力は D³ に比例 → 粒子径が2倍で抗力は8倍")

    # --- アニメーション生成 ---
    frame_dir = "/tmp/particle_frames"
    os.makedirs(frame_dir, exist_ok=True)

    fps = 15
    total_frames = 60
    steps_per_frame = sim.timesteps // total_frames   # 25 steps/frame

    # D の全体範囲を固定（カラーバーを統一）
    drag_all_max = np.max(
        3 * np.pi * sim.mu_fluid * sim.diameters
        * np.linalg.norm(np.abs(sim.v_terminal), axis=1)
    ) * 1e15 * 1.1

    print(f"\nアニメーション生成: {total_frames}フレーム × {steps_per_frame}steps …")

    def save_frame(sim, idx, drag_vmax):
        drag_fN = sim.get_drag_forces() * 1e15

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

        sc1 = axes[0].scatter(
            sim.positions[:, 0] * 1e6,
            sim.positions[:, 1] * 1e6,
            s=25 * (sim.diameters / sim.base_diameter) ** 2,
            c=drag_fN, cmap="plasma",
            vmin=0, vmax=drag_vmax,
            alpha=0.75, edgecolors="none",
        )
        axes[0].set_xlim(0, sim.domain_size * 1e6)
        axes[0].set_ylim(0, sim.domain_size * 1e6)
        axes[0].set_aspect("equal")
        axes[0].set_xlabel("X [μm]"); axes[0].set_ylabel("Y [μm]")
        axes[0].set_title(f"Stokes Drag [fN]  (proportional to D*|v|)")
        cb = plt.colorbar(sc1, ax=axes[0])
        cb.set_label("drag [fN]")

        sc2 = axes[1].scatter(
            sim.positions[:, 0] * 1e6,
            sim.positions[:, 1] * 1e6,
            s=25 * (sim.diameters / sim.base_diameter) ** 2,
            c=sim.diameters * 1e6, cmap="viridis",
            vmin=sim.base_diameter * (1 - sim.diameter_variation) * 1e6,
            vmax=sim.base_diameter * (1 + sim.diameter_variation) * 1e6,
            alpha=0.75, edgecolors="none",
        )
        axes[1].set_xlim(0, sim.domain_size * 1e6)
        axes[1].set_ylim(0, sim.domain_size * 1e6)
        axes[1].set_aspect("equal")
        axes[1].set_xlabel("X [μm]"); axes[1].set_ylabel("Y [μm]")
        axes[1].set_title("Diameter D [μm]  (larger = faster sedimentation)")
        cb2 = plt.colorbar(sc2, ax=axes[1])
        cb2.set_label("D [μm]")

        plt.suptitle(
            f"Stokes Gravitational Sedimentation  N={sim.num_particles}  "
            f"D={sim.base_diameter*1e6:.0f}um +/-{sim.diameter_variation*100:.0f}%  "
            f"t={sim.time:.3f}s  step={sim.step_count}",
            fontsize=11, fontweight="bold",
        )
        plt.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(f"{frame_dir}/frame_{idx:04d}.png", dpi=100, bbox_inches="tight")
        plt.close(fig)

    # フレーム 0: 初期状態
    save_frame(sim, 0, drag_all_max)

    for fi in range(1, total_frames + 1):
        for _ in range(steps_per_frame):
            sim.step()
        save_frame(sim, fi, drag_all_max)
        if fi % 15 == 0:
            drag_m = np.mean(sim.get_drag_forces()) * 1e15
            print(f"  frame {fi}/{total_frames}  t={sim.time:.3f}s  mean_drag={drag_m:.4f}fN")

    # ffmpeg で MP4 合成
    out_mp4 = "/tmp/particle_simulation.mp4"
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-i", f"{frame_dir}/frame_%04d.png",
        "-vcodec", "libx264", "-crf", "18",
        "-pix_fmt", "yuv420p",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        out_mp4,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        size_mb = os.path.getsize(out_mp4) / 1024 / 1024
        print(f"\n動画保存: {out_mp4}  ({size_mb:.1f} MB, {fps}fps, {total_frames/fps:.1f}s)")
    else:
        print(f"ffmpegエラー: {result.stderr[-300:]}")

    # 最終統計
    stats_f = sim.get_statistics()
    print(f"\n最終統計 (t={sim.time:.3f}s):")
    print(f"  平均速度: {stats_f['mean_velocity']*1e6:.2f} μm/s")
    print(f"  最大速度: {stats_f['max_velocity']*1e6:.2f} μm/s")
    print(f"  平均抗力: {stats_f['mean_drag_fN']:.4f} fN")
    print(f"  最大抗力: {stats_f['max_drag_fN']:.4f} fN")
    print(f"  運動エネルギー: {stats_f['kinetic_energy']:.4e} J")
    print(f"  ※ 抗力の最大/最小比 ≈ (D_max/D_min)³ = "
          f"{(stats_f['max_diameter']/stats_f['min_diameter'])**3:.1f}倍")

    # 最終フレームをプロジェクトに保存
    final_fig = sim.visualize()
    out_png = "/tmp/particle_simulation_final.png"
    final_fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(final_fig)
    print(f"\n最終状態画像: {out_png}")


if __name__ == "__main__":
    main()
