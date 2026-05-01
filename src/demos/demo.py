#!/usr/bin/env python
"""
CFD Solver ライブラリのデモンストレーション

このスクリプトは、cfd-solver ライブラリの全機能を紹介します。
"""

import matplotlib.pyplot as plt
import numpy as np

from cfd import CavityFlow, CircularPoiseuille, PlanePoiseuille


def demo_plane_poiseuille():
    """平行平板間のポアズイユ流れのデモ"""
    print("\n" + "=" * 70)
    print("デモ1: 平行平板間のポアズイユ流れ")
    print("=" * 70)
    print("""
このシミュレーションは、2枚の平行な平板の間の層流を計算します。
- 平板間距離: 10 mm
- 圧力勾配: -100 Pa/m
- 流体: 水（動粘性係数 1 mPa·s）
    """)

    # パラメータ設定
    L = 0.01  # 平板間距離 [m]
    dp_dx = -100  # 圧力勾配 [Pa/m]
    mu = 1e-3  # 動粘性係数 [Pa·s]

    # シミュレーション実行
    flow = PlanePoiseuille(L, dp_dx, mu, ny=101)

    # 解析解を計算
    u_analytical = flow.analytical_solution()

    # 数値解を計算
    u_numerical = flow.numerical_solution(tol=1e-6, max_iter=10000)

    # 結果
    u_max = flow.get_max_velocity()
    Q = flow.get_flow_rate()

    print("\n計算結果:")
    print(f"  - 最大速度（中心）: {u_max:.6f} m/s = {u_max*1000:.3f} mm/s")
    print(f"  - 流量（単位幅あたり）: {Q*1e6:.4f} mm²/s")
    print(f"  - 誤差（解析解との比較）: {np.max(np.abs(u_analytical - u_numerical)):.2e}")

    # プロット
    fig = flow.plot()
    return fig


def demo_circular_poiseuille():
    """円管内のHagen-Poiseuille流れのデモ"""
    print("\n" + "=" * 70)
    print("デモ2: 円管内のHagen-Poiseuille流れ")
    print("=" * 70)
    print("""
このシミュレーションは、円形断面の管内の層流を計算します。
- 管の半径: 5 mm
- 圧力勾配: -100 Pa/m
- 流体: 水（動粘性係数 1 mPa·s）
    """)

    # パラメータ設定
    R = 0.005  # 管の半径 [m]
    dp_dx = -100  # 圧力勾配 [Pa/m]
    mu = 1e-3  # 動粘性係数 [Pa·s]

    # シミュレーション実行
    flow = CircularPoiseuille(R, dp_dx, mu, nr=101)

    # 解析解を計算
    u_analytical = flow.analytical_solution()

    # 数値解を計算
    u_numerical = flow.numerical_solution(tol=1e-6, max_iter=10000)

    # 結果
    u_max = flow.get_max_velocity()
    Q = flow.get_flow_rate()

    print("\n計算結果:")
    print(f"  - 最大速度（中心）: {u_max:.6f} m/s = {u_max*1000:.3f} mm/s")
    print(f"  - 流量: {Q*1e9:.4f} mm³/s")
    print(f"  - 誤差（解析解との比較）: {np.max(np.abs(u_analytical - u_numerical)):.2e}")

    # プロット
    fig = flow.plot()
    return fig


def demo_cavity_flow():
    """駆動キャビティ流れのデモ"""
    print("\n" + "=" * 70)
    print("デモ3: 駆動キャビティ流れ")
    print("=" * 70)
    print("""
このシミュレーションは、正方形キャビティ内の駆動流れを計算します。
- キャビティサイズ: 1m × 1m
- 上蓋速度: 1 m/s
- 流体密度: 1 kg/m³
- 動粘性係数: 0.01 Pa·s
- レイノルズ数: 100
    """)

    # パラメータ設定
    L = 1.0
    rho = 1.0
    mu = 0.01

    # シミュレーション実行
    cavity = CavityFlow(L, rho, mu, nx=65, ny=65)

    # レイノルズ数を表示
    Re = cavity.get_Reynolds_number()
    print("\nシミュレーション情報:")
    print(f"  - レイノルズ数: Re = {Re:.2f}")
    print(f"  - グリッド解像度: {cavity.nx} × {cavity.ny}")
    print(f"  - 上蓋速度: {cavity.U_lid} m/s")

    print("\n計算中...")
    cavity.solve_steady_state(max_iterations=300, verbose=False)
    print("計算完了！")

    # プロット
    figs = []
    print("\n可視化生成中...")
    figs.append(cavity.plot_velocity_field())
    figs.append(cavity.plot_streamlines())
    figs.append(cavity.plot_pressure_field())

    return figs


def print_capabilities():
    """ライブラリの機能を表示"""
    print("\n" + "=" * 70)
    print("CFD Solver ライブラリの機能一覧")
    print("=" * 70)

    print("""
【ポアズイユ流れ (Poiseuille Flow)】

1. 平行平板間の流れ (PlanePoiseuille)
   - 2枚の平行な平板間の層流シミュレーション
   - 支配方程式: d²u/dy² = (1/μ)·dp/dx
   - 解析解: パラボラ分布
   - 計算方法: 有限差分法 + ガウス・ザイデル反復法
   
   できること:
   ✓ 速度分布の計算（解析解・数値解）
   ✓ 最大速度の計算
   ✓ 流量（単位幅あたり）の計算
   ✓ 解析解と数値解の比較
   ✓ 可視化（速度プロファイル）

2. 円管内の流れ (CircularPoiseuille)
   - 円形断面の管内の層流シミュレーション（Hagen-Poiseuille流れ）
   - 支配方程式: (1/r)·d/dr(r·du/dr) = (1/μ)·dp/dx
   - 解析解: パラボラ分布（円形）
   - 計算方法: 有限差分法（円筒座標系）+ ガウス・ザイデル反復法
   
   できること:
   ✓ 速度分布の計算（解析解・数値解）
   ✓ 最大速度の計算
   ✓ Hagen-Poiseuille流量公式による流量計算
   ✓ 解析解と数値解の比較
   ✓ 1D・2D可視化（速度プロファイル・断面図）

【キャビティ流れ (Cavity Flow)】

3. 駆動キャビティ流れ (CavityFlow)
   - 正方形キャビティ内の駆動層流シミュレーション
   - 支配方程式: 2次元非圧縮性ナビエ・ストークス方程式
   - 計算方法: 陽的オイラー法 + SOR法（圧力求解）
   
   できること:
   ✓ 非圧縮性流れの数値シミュレーション
   ✓ ナビエ・ストークス方程式の求解
   ✓ 定常状態までの時間発展計算
   ✓ レイノルズ数の計算
   ✓ 速度場の可視化（速度ベクトル・スカラー場）
   ✓ 流線の可視化
   ✓ 圧力場の可視化
   ✓ 渦度計算
   ✓ 発散（連続方程式残差）の計算

【共通機能】

✓ 解析的な結果と数値的な結果の比較
✓ 物理量の計算（速度、流量、圧力など）
✓ matplotlib によるプロット・可視化
✓ NumPy による高速計算
✓ パラメータの柔軟な設定

【応用例】

1. 教育用途
   - 流体力学の基礎概念の理解
   - 数値計算方法の学習
   - 有限差分法の実装例

2. 研究用途
   - パラメータスタディ
   - 異なる条件での流れ特性の検証
   - 数値解析法の精度検証

3. 設計支援
   - 流路設計の初期検討
   - フロー分析
   - 性能予測
    """)


def main():
    """メイン実行"""
    print("\n" + "=" * 70)
    print("CFD Solver - 計算流体力学ライブラリ デモンストレーション")
    print("=" * 70)

    # ライブラリの機能を表示
    print_capabilities()

    # デモンストレーション実行
    figs = []

    # デモ1: 平行平板間
    try:
        fig1 = demo_plane_poiseuille()
        figs.append(fig1)
    except Exception as e:
        print(f"デモ1 エラー: {e}")

    # デモ2: 円管内
    try:
        fig2 = demo_circular_poiseuille()
        figs.append(fig2)
    except Exception as e:
        print(f"デモ2 エラー: {e}")

    # デモ3: キャビティ流れ
    try:
        figs3 = demo_cavity_flow()
        figs.extend(figs3)
    except Exception as e:
        print(f"デモ3 エラー: {e}")

    # プロットを表示
    print("\n" + "=" * 70)
    print("すべてのデモ実行完了！")
    print("=" * 70)
    print(f"生成された図：{len(figs)} 個")

    # 表示（Jupyter環境の場合）
    try:
        plt.show()
    except Exception:
        print("（プロット表示はスキップされました）")


if __name__ == "__main__":
    main()
