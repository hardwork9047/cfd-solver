#!/usr/bin/env python
"""
キャビティ流れの診断スクリプト

現在のアルゴリズムの問題点を診断し、改善案を提案します。
"""

import numpy as np
import matplotlib.pyplot as plt
from cfd import CavityFlow


def diagnose_cavity_flow():
    """キャビティ流れの問題診断"""
    print("\n" + "="*70)
    print("キャビティ流れ診断")
    print("="*70)
    
    # 小さなグリッドでテスト
    print("\n[診断1] グリッド解像度の影響")
    print("-" * 70)
    
    for nx in [17, 33]:
        ny = nx
        cavity = CavityFlow(L=1.0, rho=1.0, mu=0.01, nx=nx, ny=ny)
        Re = cavity.get_Reynolds_number()
        
        print(f"\nグリッド: {nx}×{ny}, Re={Re:.1f}")
        print(f"  dt={cavity.dt:.2e}, dx={cavity.dx:.4f}, dy={cavity.dy:.4f}")
        print(f"  nu={cavity.nu:.4f} (kinematic viscosity)")
        
        # 数ステップ実行
        print("  計算中...", end="", flush=True)
        try:
            for i in range(20):
                cavity.step()
            print(" OK")
            
            # 統計
            u_max = np.max(np.abs(cavity.u))
            v_max = np.max(np.abs(cavity.v))
            p_max = np.max(np.abs(cavity.p))
            
            print(f"  u_max={u_max:.4f}, v_max={v_max:.4f}, p_max={p_max:.2e}")
            
            # NaN チェック
            if np.any(np.isnan(cavity.u)) or np.any(np.isnan(cavity.v)):
                print("  ⚠️  NaN値が発生！")
            else:
                print("  ✓ NaN なし")
                
        except Exception as e:
            print(f" ERROR: {e}")
    
    # より詳細な分析
    print("\n[診断2] 時間発展の安定性")
    print("-" * 70)
    
    cavity = CavityFlow(L=1.0, rho=1.0, mu=0.01, nx=33, ny=33)
    
    u_history = []
    v_history = []
    p_history = []
    
    print("計算中...", end="", flush=True)
    for step in range(50):
        cavity.step()
        u_history.append(np.max(np.abs(cavity.u)))
        v_history.append(np.max(np.abs(cavity.v)))
        p_history.append(np.max(np.abs(cavity.p)))
        
        if step % 10 == 0:
            print(f".", end="", flush=True)
    
    print(" 完了")
    
    # プロット
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    
    ax1.semilogy(u_history, 'b-o', linewidth=2, markersize=4)
    ax1.set_xlabel('Step')
    ax1.set_ylabel('max|u| [m/s]')
    ax1.set_title('x方向速度の時間変化')
    ax1.grid(True, alpha=0.3)
    
    ax2.semilogy(v_history, 'r-o', linewidth=2, markersize=4)
    ax2.set_xlabel('Step')
    ax2.set_ylabel('max|v| [m/s]')
    ax2.set_title('y方向速度の時間変化')
    ax2.grid(True, alpha=0.3)
    
    ax3.semilogy(p_history, 'g-o', linewidth=2, markersize=4)
    ax3.set_xlabel('Step')
    ax3.set_ylabel('max|p| [Pa]')
    ax3.set_title('圧力の時間変化')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/tmp/cavity_stability.png', dpi=100, bbox_inches='tight')
    print(f"\n保存: /tmp/cavity_stability.png")
    
    # 最終的な速度場を表示
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # u分布
    im1 = ax1.contourf(cavity.X, cavity.Y, cavity.u, levels=20, cmap='RdBu_r')
    ax1.set_title('u速度分布')
    ax1.set_aspect('equal')
    plt.colorbar(im1, ax=ax1)
    
    # v分布
    im2 = ax2.contourf(cavity.X, cavity.Y, cavity.v, levels=20, cmap='RdBu_r')
    ax2.set_title('v速度分布')
    ax2.set_aspect('equal')
    plt.colorbar(im2, ax=ax2)
    
    plt.tight_layout()
    plt.savefig('/tmp/cavity_velocity.png', dpi=100, bbox_inches='tight')
    print(f"保存: /tmp/cavity_velocity.png")
    
    # 問題点の分析
    print("\n[診断結果] 検出された問題")
    print("-" * 70)
    
    if u_history[-1] < 0.01:
        print("❌ 速度が小さすぎる（収束がない、または発散している）")
    elif u_history[-1] > 2:
        print("❌ 速度が大きすぎる（数値不安定性の可能性）")
    else:
        print("✓ 速度の大きさは適切な範囲内")
    
    # 増減傾向の確認
    u_growth = (u_history[-1] - u_history[0]) / (u_history[0] + 1e-10)
    print(f"\n速度成長率: {u_growth:.2%}")
    
    if u_growth > 0:
        print("✓ 速度が増加（初期状態から流れが発展）")
    else:
        print("❌ 速度が減少（数値粘性が強すぎる？）")


if __name__ == "__main__":
    diagnose_cavity_flow()
