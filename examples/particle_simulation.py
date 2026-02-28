#!/usr/bin/env python
"""
粒子シミュレーション：粒子径を乱数で変化させる

パラメータ：
- 粒子数: 100個
- 粒子径: 基本値 ±15% の範囲で乱数変化
- タイムステップ数: 5倍（300ステップ）
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

# ロギング設定
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


class ParticleSimulation:
    """乱数粒子径を持つ粒子シミュレーション"""

    def __init__(
        self,
        num_particles=100,
        base_diameter=1.0e-6,  # 基本粒子径 [m]
        diameter_variation=0.15,  # 粒子径変動 ±15%
        domain_size=1.0e-4,  # シミュレーション領域 [m]
        timesteps=300,  # タイムステップ数 (5倍)
        dt=1.0e-5,  # 時間ステップ [s]
    ):
        """
        初期化

        Parameters:
        -----------
        num_particles : int
            粒子数
        base_diameter : float
            基本粒子径 [m]
        diameter_variation : float
            粒子径の変動幅（±の割合）
        domain_size : float
            シミュレーション領域サイズ [m]
        timesteps : int
            タイムステップ数
        dt : float
            時間ステップサイズ [s]
        """
        self.num_particles = num_particles
        self.base_diameter = base_diameter
        self.diameter_variation = diameter_variation
        self.domain_size = domain_size
        self.timesteps = timesteps
        self.dt = dt

        # 粒子座標 [m]
        self.positions = np.random.uniform(
            0, domain_size, size=(num_particles, 2)
        )

        # 粒子速度 [m/s]
        self.velocities = np.random.normal(
            0, 0.01, size=(num_particles, 2)
        )

        # 粒子径を乱数で生成（±15% の範囲）
        # 均一分布: [base_diameter * (1 - variation), base_diameter * (1 + variation)]
        min_diameter = base_diameter * (1 - diameter_variation)
        max_diameter = base_diameter * (1 + diameter_variation)
        self.diameters = np.random.uniform(
            min_diameter, max_diameter, size=num_particles
        )

        # 粒子の質量（密度を一定と�定）
        density = 1000  # kg/m³ (水の密度)
        volume = (4 / 3) * np.pi * (self.diameters / 2) ** 3
        self.masses = density * volume

        # 加速度（重力など）
        self.accelerations = np.zeros((num_particles, 2))
        self.accelerations[:, 1] = -9.8  # 重力加速度 [m/s²]

        self.time = 0.0
        self.step_count = 0

        logger.info(f"粒子シミュレーション初期化:")
        logger.info(f"  粒子数: {num_particles}")
        logger.info(f"  基本粒子径: {base_diameter*1e6:.2f} μm")
        logger.info(f"  粒子径範囲: {min_diameter*1e6:.2f} - {max_diameter*1e6:.2f} μm")
        logger.info(f"  タイムステップ数: {timesteps}")
        logger.info(f"  時間ステップサイズ: {dt:.2e} s")

    def step(self):
        """1ステップ計算"""
        # 速度の更新 (v = v + a*dt)
        self.velocities += self.accelerations * self.dt

        # 位置の更新 (x = x + v*dt)
        self.positions += self.velocities * self.dt

        # 壁との衝突判定（簡単な反射）
        for i in range(self.num_particles):
            # 下の壁
            if self.positions[i, 1] < 0:
                self.positions[i, 1] = 0
                self.velocities[i, 1] *= -0.8  # 反発係数 0.8

            # 上の壁
            if self.positions[i, 1] > self.domain_size:
                self.positions[i, 1] = self.domain_size
                self.velocities[i, 1] *= -0.8

            # 左の壁
            if self.positions[i, 0] < 0:
                self.positions[i, 0] = 0
                self.velocities[i, 0] *= -0.8

            # 右の壁
            if self.positions[i, 0] > self.domain_size:
                self.positions[i, 0] = self.domain_size
                self.velocities[i, 0] *= -0.8

        self.time += self.dt
        self.step_count += 1

    def run(self):
        """シミュレーションを実行"""
        logger.info("シミュレーション実行中...")

        for step in range(self.timesteps):
            self.step()
            if (step + 1) % 50 == 0:
                logger.info(f"  ステップ {step + 1}/{self.timesteps} - 時刻: {self.time:.4e} s")

        logger.info("シミュレーション完了！")

    def visualize(self):
        """粒子の現在状態を可視化"""
        fig, ax = plt.subplots(figsize=(10, 10))

        # 粒子を円として描画
        scatter = ax.scatter(
            self.positions[:, 0],
            self.positions[:, 1],
            s=50 * (self.diameters / self.base_diameter) ** 2,  # サイズを粒子径に対応
            c=self.diameters * 1e6,  # 色を粒子径で表示
            cmap="viridis",
            alpha=0.6,
            edgecolors="black",
            linewidth=0.5,
        )

        ax.set_xlim(0, self.domain_size)
        ax.set_ylim(0, self.domain_size)
        ax.set_aspect("equal")
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")
        ax.set_title(
            f"粒子シミュレーション\n時刻: {self.time:.4e} s, ステップ: {self.step_count}"
        )

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("粒子径 [μm]")

        plt.tight_layout()
        return fig

    def get_statistics(self):
        """粒子の統計情報を取得"""
        v_mag = np.linalg.norm(self.velocities, axis=1)

        stats = {
            "mean_diameter": np.mean(self.diameters) * 1e6,  # μm
            "std_diameter": np.std(self.diameters) * 1e6,   # μm
            "min_diameter": np.min(self.diameters) * 1e6,   # μm
            "max_diameter": np.max(self.diameters) * 1e6,   # μm
            "mean_velocity": np.mean(v_mag),               # m/s
            "max_velocity": np.max(v_mag),                 # m/s
            "kinetic_energy": 0.5 * np.sum(self.masses * v_mag ** 2),  # J
        }

        return stats


def main():
    """メイン実行"""
    print("\n" + "=" * 70)
    print("粒子シミュレーション：乱数粒子径")
    print("=" * 70)

    # パラメータ
    num_particles = 100
    base_diameter = 1.0e-6  # 1 μm
    diameter_variation = 0.15  # ±15%
    domain_size = 1.0e-4  # 100 μm
    timesteps = 300  # 5倍のタイムステップ
    dt = 1.0e-5  # 10 μs

    # シミュレーション実行
    sim = ParticleSimulation(
        num_particles=num_particles,
        base_diameter=base_diameter,
        diameter_variation=diameter_variation,
        domain_size=domain_size,
        timesteps=timesteps,
        dt=dt,
    )

    print(f"\n初期粒子径統計:")
    stats_init = sim.get_statistics()
    print(f"  平均粒子径: {stats_init['mean_diameter']:.3f} μm")
    print(f"  標準偏差: {stats_init['std_diameter']:.3f} μm")
    print(f"  粒子径範囲: {stats_init['min_diameter']:.3f} - {stats_init['max_diameter']:.3f} μm")

    # シミュレーション実行
    sim.run()

    # 統計情報を出力
    print(f"\n最終統計:")
    stats_final = sim.get_statistics()
    print(f"  平均粒子径: {stats_final['mean_diameter']:.3f} μm")
    print(f"  標準偏差: {stats_final['std_diameter']:.3f} μm")
    print(f"  平均速度: {stats_final['mean_velocity']:.4f} m/s")
    print(f"  最大速度: {stats_final['max_velocity']:.4f} m/s")
    print(f"  運動エネルギー: {stats_final['kinetic_energy']:.4e} J")

    # 可視化
    print("\n可視化生成中...")
    fig = sim.visualize()
    plt.savefig("/tmp/particle_simulation_final.png", dpi=150, bbox_inches="tight")
    print("保存: /tmp/particle_simulation_final.png")

    plt.show()


if __name__ == "__main__":
    main()
