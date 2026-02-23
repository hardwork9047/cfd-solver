#!/usr/bin/env python
"""
改善されたキャビティ流れのデモンストレーション

低レイノルズ数から始めて、徐々に複雑な流れへ
"""

import matplotlib.pyplot as plt
import numpy as np

from cfd import CavityFlow


def demo_low_re_cavity():
    """低レイノルズ数でのキャビティ流れ"""
    print("\n" + "=" * 70)
    print("デモ: 低レイノルズ数（Re=10）でのキャビティ流れ")
    print("=" * 70)

    # Re = 10 の条件を設定
    # Re = ρ·U·L / μ = 10
    mu = 1.0 / 10  # ρ=1, U=1, L=1 の場合

    cavity = CavityFlow(L=1.0, rho=1.0, mu=mu, nx=49, ny=49)
    Re = cavity.get_Reynolds_number()

    print(f"レイノルズ数: Re = {Re:.1f}")
    print(f"グリッド: {cavity.nx}×{cavity.ny}")
    print(f"動粘性係数: μ = {mu:.4f} Pa·s")
    print(f"初期 dt: {cavity.dt:.2e}")
    print()

    print("定常状態へ計算中...", end="", flush=True)
    cavity.solve_steady_state(max_iterations=500, verbose=False)
    print(" 完了")
    print()

    print(f"計算時間: {cavity.time:.4f}s")
    print("最終状態:")
    print(f"  - u_max: {np.max(np.abs(cavity.u)):.4f} m/s")
    print(f"  - v_max: {np.max(np.abs(cavity.v)):.4f} m/s")
    print(f"  - p_max: {np.max(np.abs(cavity.p)):.2e} Pa")

    # 可視化
    print("\n可視化生成中...")

    fig1 = cavity.plot_velocity_field()
    fig1.savefig("/tmp/cavity_low_re_velocity.png", dpi=100, bbox_inches="tight")
    print("保存: /tmp/cavity_low_re_velocity.png")

    fig2 = cavity.plot_streamlines()
    fig2.savefig("/tmp/cavity_low_re_streamlines.png", dpi=100, bbox_inches="tight")
    print("保存: /tmp/cavity_low_re_streamlines.png")

    fig3 = cavity.plot_pressure_field()
    fig3.savefig("/tmp/cavity_low_re_pressure.png", dpi=100, bbox_inches="tight")
    print("保存: /tmp/cavity_low_re_pressure.png")

    return cavity


def demo_mid_re_cavity():
    """中程度のレイノルズ数（Re=100）でのキャビティ流れ"""
    print("\n" + "=" * 70)
    print("デモ: 中程度のレイノルズ数（Re=100）でのキャビティ流れ")
    print("=" * 70)

    # Re = 100 の条件を設定
    mu = 1.0 / 100  # ρ=1, U=1, L=1 の場合

    cavity = CavityFlow(L=1.0, rho=1.0, mu=mu, nx=65, ny=65)
    Re = cavity.get_Reynolds_number()

    print(f"レイノルズ数: Re = {Re:.1f}")
    print(f"グリッド: {cavity.nx}×{cavity.ny}")
    print(f"動粘性係数: μ = {mu:.6f} Pa·s")
    print(f"初期 dt: {cavity.dt:.2e}")
    print()

    print("定常状態へ計算中...", end="", flush=True)
    cavity.solve_steady_state(max_iterations=1000, verbose=False)
    print(" 完了")
    print()

    print(f"計算時間: {cavity.time:.4f}s")
    print("最終状態:")
    print(f"  - u_max: {np.max(np.abs(cavity.u)):.4f} m/s")
    print(f"  - v_max: {np.max(np.abs(cavity.v)):.4f} m/s")
    print(f"  - p_max: {np.max(np.abs(cavity.p)):.2e} Pa")

    # 可視化
    print("\n可視化生成中...")

    fig1 = cavity.plot_velocity_field()
    fig1.savefig("/tmp/cavity_mid_re_velocity.png", dpi=100, bbox_inches="tight")
    print("保存: /tmp/cavity_mid_re_velocity.png")

    fig2 = cavity.plot_streamlines()
    fig2.savefig("/tmp/cavity_mid_re_streamlines.png", dpi=100, bbox_inches="tight")
    print("保存: /tmp/cavity_mid_re_streamlines.png")

    return cavity


def plot_comparison(cavity_low, cavity_mid):
    """2つのシミュレーション結果を比較"""
    print("\n" + "=" * 70)
    print("レイノルズ数の影響比較")
    print("=" * 70)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Re=10 の速度
    im1 = axes[0, 0].contourf(
        cavity_low.X,
        cavity_low.Y,
        np.sqrt(cavity_low.u**2 + cavity_low.v**2),
        levels=20,
        cmap="viridis",
    )
    axes[0, 0].set_title("Velocity Magnitude (Re=10)", fontsize=12, fontweight="bold")
    axes[0, 0].set_aspect("equal")
    plt.colorbar(im1, ax=axes[0, 0])

    # Re=100 の速度
    im2 = axes[0, 1].contourf(
        cavity_mid.X,
        cavity_mid.Y,
        np.sqrt(cavity_mid.u**2 + cavity_mid.v**2),
        levels=20,
        cmap="viridis",
    )
    axes[0, 1].set_title("Velocity Magnitude (Re=100)", fontsize=12, fontweight="bold")
    axes[0, 1].set_aspect("equal")
    plt.colorbar(im2, ax=axes[0, 1])

    # Re=10 の流線
    axes[1, 0].streamplot(
        cavity_low.X, cavity_low.Y, cavity_low.u, cavity_low.v, color="steelblue", density=1.5
    )
    axes[1, 0].set_title("Streamlines (Re=10)", fontsize=12, fontweight="bold")
    axes[1, 0].set_aspect("equal")

    # Re=100 の流線
    axes[1, 1].streamplot(
        cavity_mid.X, cavity_mid.Y, cavity_mid.u, cavity_mid.v, color="steelblue", density=1.5
    )
    axes[1, 1].set_title("Streamlines (Re=100)", fontsize=12, fontweight="bold")
    axes[1, 1].set_aspect("equal")

    plt.tight_layout()
    plt.savefig("/tmp/cavity_comparison.png", dpi=100, bbox_inches="tight")
    print("\n保存: /tmp/cavity_comparison.png")


def main():
    """メイン実行"""
    print("\n" + "=" * 70)
    print("改善版キャビティ流れシミュレーション - デモンストレーション")
    print("=" * 70)

    # デモ1: 低レイノルズ数
    try:
        cavity_low = demo_low_re_cavity()
    except Exception as e:
        print(f"❌ デモ1 エラー: {e}")
        import traceback

        traceback.print_exc()
        cavity_low = None

    # デモ2: 中程度のレイノルズ数
    try:
        cavity_mid = demo_mid_re_cavity()
    except Exception as e:
        print(f"❌ デモ2 エラー: {e}")
        import traceback

        traceback.print_exc()
        cavity_mid = None

    # 比較
    if cavity_low is not None and cavity_mid is not None:
        try:
            plot_comparison(cavity_low, cavity_mid)
        except Exception as e:
            print(f"❌ 比較プロット エラー: {e}")

    print("\n" + "=" * 70)
    print("デモンストレーション完了！")
    print("=" * 70)


if __name__ == "__main__":
    main()
