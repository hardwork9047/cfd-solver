"""
LBM-DEM 連成シミュレーション 実行スクリプト
==============================================
- 流体: D2Q9 LBM (Poiseuille チャネル流)
- 粒子: DEM (Hertz接触 + Stokes抗力)
- 動画: MP4 (ffmpeg) を出力
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np

# パッケージを src/ から読み込む
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from cfd_dem_lbm import LBMDEMSolver

# ---------------------------------------------------------------------------
# コマンドライン引数
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="LBM-DEM coupled simulation runner")
parser.add_argument("--cylinder", action="store_true", help="Add a fixed solid cylinder")
parser.add_argument("--cyl-x", type=float, default=None,
                    help="Cylinder center x [lattice] (default: NX/4)")
parser.add_argument("--cyl-y", type=float, default=None,
                    help="Cylinder center y [lattice] (default: NY/2)")
parser.add_argument("--cyl-r", type=float, default=6.0,
                    help="Cylinder radius [lattice] (default: 6.0)")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# キャッシュ付きサブクラス: _macroscopic() を各ステップ1回だけ計算
# ---------------------------------------------------------------------------

class FastLBMDEM(LBMDEMSolver):
    """マクロ量をキャッシュして高速化したサブクラス。"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cached_rho = None
        self._cached_ux = None
        self._cached_uy = None
        self._cache_step = -1

    def _refresh_cache(self):
        self._cached_rho, self._cached_ux, self._cached_uy = super()._macroscopic()
        self._cache_step = self.step_count

    def _macroscopic(self):
        if self._cache_step != self.step_count:
            self._refresh_cache()
        return self._cached_rho, self._cached_ux, self._cached_uy

    def advance(self, n_steps: int = 1):
        for _ in range(n_steps):
            # LBM ステップ前にキャッシュを更新
            self._refresh_cache()
            # 親クラスの advance を1ステップ呼び出す
            super().advance(1)


# ---------------------------------------------------------------------------
# シミュレーション設定
# ---------------------------------------------------------------------------

NX              = 180    # 格子幅
NY              = 70     # 格子高さ
RE              = 100.0  # Reynolds 数
U_MAX           = 0.05   # 最大流速 (格子単位)
N_PARTICLES     = 40     # 粒子数
RADIUS          = 3.0    # 粒子平均半径 (格子単位)
RADIUS_VARIATION = 0.15  # 粒子径バリエーション ±15%
DENSITY_RATIO   = 2.0    # ρ_p / ρ_f
GRAVITY         = 3e-5   # 重力加速度 (格子単位)
TOTAL_STEPS  = 20000  # 総 LBM ステップ数
SNAPSHOT_EVERY = 40   # 何ステップごとにスナップショット取得
FPS          = 24     # 動画フレームレート
OUT_DIR      = Path("/tmp/lbm_dem_results")

OUT_DIR.mkdir(exist_ok=True)

# 円柱設定 (--cylinder が指定された場合のみ有効)
if args.cylinder:
    CYL_X = args.cyl_x if args.cyl_x is not None else NX / 4
    CYL_Y = args.cyl_y if args.cyl_y is not None else NY / 2
    CYL_R = args.cyl_r
    CYLINDER = (CYL_X, CYL_Y, CYL_R)
else:
    CYLINDER = None

# ---------------------------------------------------------------------------
# 初期化
# ---------------------------------------------------------------------------

print("=" * 60)
print("  LBM-DEM 連成シミュレーション")
print("=" * 60)

sim = FastLBMDEM(
    nx=NX, ny=NY,
    Re=RE, u_max=U_MAX,
    n_particles=N_PARTICLES,
    particle_radius=RADIUS,
    radius_variation=RADIUS_VARIATION,
    density_ratio=DENSITY_RATIO,
    gravity=GRAVITY,
    dem_substeps=4,
    seed=42,
    cylinder=CYLINDER,
)

# ---------------------------------------------------------------------------
# ウォームアップ (流れを発達させる)
# ---------------------------------------------------------------------------

WARMUP = 1000
print(f"\n[1/3] ウォームアップ {WARMUP} ステップ ...")
t0 = time.perf_counter()
sim.advance(WARMUP)
print(f"      完了 ({time.perf_counter()-t0:.1f} s)")

# ---------------------------------------------------------------------------
# スナップショットを取りながらシミュレーション
# ---------------------------------------------------------------------------

n_frames = TOTAL_STEPS // SNAPSHOT_EVERY
print(f"\n[2/3] メインシミュレーション {TOTAL_STEPS} ステップ "
      f"(スナップショット {n_frames} 枚) ...")

snapshots = []   # list of (rho, ux, uy, pos_copy, vel_copy, step)

t1 = time.perf_counter()
for frame_idx in range(n_frames):
    sim.advance(SNAPSHOT_EVERY)

    rho, ux, uy = sim.get_fields()

    # 各粒子に加わる合力の大きさ (重力 + 抗力 + 接触力)
    total_force_mags = np.linalg.norm(sim.forces_p, axis=1)

    snapshots.append({
        "step": sim.step_count,
        "ux":   ux.copy(),
        "uy":   uy.copy(),
        "speed": np.sqrt(ux**2 + uy**2).copy(),
        "pos":  sim.pos.copy(),
        "vel":  sim.vel.copy(),
        "total_force": total_force_mags.copy(),
    })

    elapsed = time.perf_counter() - t1
    frac = (frame_idx + 1) / n_frames
    eta = elapsed / frac - elapsed if frac > 0 else 0
    speed_max = float(np.sqrt(ux**2 + uy**2).max())
    p_ke = 0.5 * sim.mass_p * float(np.sum(sim.vel**2))
    f_mean = total_force_mags.mean()
    print(f"  frame {frame_idx+1:>3}/{n_frames}  step={sim.step_count:>6,}"
          f"  |u|_max={speed_max:.5f}  KE_p={p_ke:.3e}  |F|_mean={f_mean:.3e}  ETA {eta:.0f}s")

total_time = time.perf_counter() - t1
print(f"      完了 ({total_time:.1f} s, {TOTAL_STEPS/total_time:.0f} steps/s)")

# ---------------------------------------------------------------------------
# 統計出力
# ---------------------------------------------------------------------------

print("\n--- 最終統計 ---")
last = snapshots[-1]
p_vel_mag = np.linalg.norm(last["vel"], axis=1)
print(f"  粒子数          : {sim.n_p}")
print(f"  流体最大速度    : {last['speed'].max():.5f} (格子単位)")
print(f"  粒子平均速度    : {p_vel_mag.mean():.4e} (格子単位)")
print(f"  粒子最大速度    : {p_vel_mag.max():.4e} (格子単位)")
print(f"  粒子KE          : {0.5*sim.mass_p*np.sum(last['vel']**2):.4e}")
print(f"  粒子Y重心       : {last['pos'][:,1].mean():.2f} / {NY} (格子単位)")
print(f"  粒子Y重心 (正規): {last['pos'][:,1].mean()/NY:.3f}")

# 静止画を保存
fig_stat, axes = plt.subplots(1, 3, figsize=(16, 5))
fig_stat.suptitle(
    f"LBM-DEM  Re={RE:.0f}  step={sim.step_count:,}  {sim.n_p} particles  ({NX}×{NY})",
    fontsize=12,
)
snap = snapshots[-1]
x = np.arange(NX)
y = np.arange(NY)

def _add_cylinder_patch(ax, color="cyan", alpha=0.6, zorder=4):
    """静止画用: 円柱パッチを追加する。"""
    if CYLINDER is not None:
        ax.add_patch(plt.Circle(
            (CYLINDER[0], CYLINDER[1]), CYLINDER[2],
            color=color, alpha=alpha, zorder=zorder,
        ))

ax = axes[0]
im = ax.imshow(snap["speed"].T, origin="lower", cmap="inferno",
               extent=[0, NX, 0, NY], aspect="auto")
for i in range(sim.n_p):
    c = plt.Circle((snap["pos"][i,0], snap["pos"][i,1]), RADIUS,
                   color="white", lw=0.8, fill=False)
    ax.add_patch(c)
_add_cylinder_patch(ax)
ax.set_title("速度大きさ |u|")
ax.set_xlabel("x [格子]"); ax.set_ylabel("y [格子]")
fig_stat.colorbar(im, ax=ax, shrink=0.7)

ax = axes[1]
lw_arr = 1.5 * snap["speed"].T / (snap["speed"].T.max() + 1e-12)
ax.streamplot(x, y, snap["ux"].T, snap["uy"].T,
              color=snap["speed"].T, cmap="cool",
              linewidth=lw_arr, density=1.2, arrowsize=0.8)
for i in range(sim.n_p):
    c = plt.Circle((snap["pos"][i,0], snap["pos"][i,1]), RADIUS,
                   color="white", lw=0.8, fill=False)
    ax.add_patch(c)
_add_cylinder_patch(ax)
ax.set_xlim(0, NX); ax.set_ylim(0, NY)
ax.set_title("流線"); ax.set_xlabel("x [格子]")

ax = axes[2]
sc = ax.scatter(snap["pos"][:,0], snap["pos"][:,1],
                c=snap["total_force"], cmap="plasma", s=(RADIUS*4)**2,
                edgecolors="k", lw=0.5, zorder=5)
ax.imshow(snap["speed"].T, origin="lower", cmap="Blues",
          extent=[0, NX, 0, NY], aspect="auto", alpha=0.5)
_add_cylinder_patch(ax)
fig_stat.colorbar(sc, ax=ax, label="合力 |F_total| [格子単位]", shrink=0.7)
ax.set_xlim(0, NX); ax.set_ylim(0, NY)
ax.set_title("粒子位置 (合力でカラー)"); ax.set_xlabel("x [格子]")

plt.tight_layout()
static_path = OUT_DIR / "lbm_dem_final.png"
fig_stat.savefig(static_path, dpi=150, bbox_inches="tight")
print(f"\n静止画保存: {static_path}")
plt.close(fig_stat)

# ---------------------------------------------------------------------------
# 動画作成
# ---------------------------------------------------------------------------

print(f"\n[3/3] 動画作成 ({n_frames} フレーム, {FPS} fps) ...")

fig_anim, ax_anim = plt.subplots(figsize=(10, 4))
fig_anim.patch.set_facecolor("#0a0a0a")
ax_anim.set_facecolor("#0a0a0a")

snap0 = snapshots[0]
speed_global_max = max(s["speed"].max() for s in snapshots)

# 合力の全フレームにわたるグローバルmin/max（カラースケール固定）
force_global_min = min(s["total_force"].min() for s in snapshots)
force_global_max = max(s["total_force"].max() for s in snapshots)
force_cmap = plt.cm.plasma
force_norm = plt.Normalize(vmin=force_global_min, vmax=force_global_max)

im_fluid = ax_anim.imshow(
    snap0["speed"].T,
    origin="lower", cmap="inferno",
    extent=[0, NX, 0, NY], aspect="auto",
    vmin=0, vmax=speed_global_max,
    animated=True,
)
cbar_fluid = fig_anim.colorbar(im_fluid, ax=ax_anim, shrink=0.75, pad=0.01)
cbar_fluid.set_label("|u| [格子単位]", color="white")
cbar_fluid.ax.yaxis.set_tick_params(color="white")
plt.setp(cbar_fluid.ax.yaxis.get_ticklabels(), color="white")

# 合力カラーバー（ScalarMappable で追加）
sm_force = plt.cm.ScalarMappable(cmap=force_cmap, norm=force_norm)
sm_force.set_array([])
cbar_force = fig_anim.colorbar(sm_force, ax=ax_anim, shrink=0.75, pad=0.12)
cbar_force.set_label("|F_total| [格子単位]", color="white")
cbar_force.ax.yaxis.set_tick_params(color="white")
plt.setp(cbar_force.ax.yaxis.get_ticklabels(), color="white")

# 固定円柱パッチ (アニメーション中は動かないので animated=False)
if CYLINDER is not None:
    ax_anim.add_patch(mpatches.Circle(
        (CYLINDER[0], CYLINDER[1]), CYLINDER[2],
        linewidth=1.5, edgecolor="cyan", facecolor="cyan", alpha=0.55, zorder=4,
    ))

# 粒子の円パッチ (合力に応じた初期色、per-particle半径でサイズ)
particle_circles = []
for i in range(sim.n_p):
    rgba = force_cmap(force_norm(snap0["total_force"][i]))
    c = mpatches.Circle(
        (snap0["pos"][i, 0], snap0["pos"][i, 1]),
        sim.radii[i],
        linewidth=1.2,
        edgecolor="white",
        facecolor=rgba,
        alpha=0.85,
        animated=True,
    )
    ax_anim.add_patch(c)
    particle_circles.append(c)

title_txt = ax_anim.set_title(
    f"LBM-DEM  step={snap0['step']:,}", color="white", fontsize=11
)
ax_anim.set_xlabel("x [格子]", color="white")
ax_anim.set_ylabel("y [格子]", color="white")
ax_anim.tick_params(colors="white")
for spine in ax_anim.spines.values():
    spine.set_edgecolor("white")

plt.tight_layout()


def update(frame_idx: int):
    snap = snapshots[frame_idx]
    im_fluid.set_data(snap["speed"].T)
    for i, c in enumerate(particle_circles):
        c.center = (snap["pos"][i, 0], snap["pos"][i, 1])
        # 合力の大きさに応じて色を更新
        rgba = force_cmap(force_norm(snap["total_force"][i]))
        c.set_facecolor(rgba)
    p_ke = 0.5 * sim.mass_p * float(np.sum(snap["vel"] ** 2))
    force_mean = snap["total_force"].mean()
    title_txt.set_text(
        f"LBM-DEM  Re={RE:.0f}  step={snap['step']:,}"
        f"  KE_p={p_ke:.3e}  |F_total|_mean={force_mean:.3e}"
    )
    artists = [im_fluid, title_txt] + particle_circles
    return artists


ani = animation.FuncAnimation(
    fig_anim, update,
    frames=n_frames,
    interval=1000 / FPS,
    blit=True,
)

video_path = OUT_DIR / "lbm_dem_simulation.mp4"
writer = animation.FFMpegWriter(fps=FPS, bitrate=2000,
                                 extra_args=["-pix_fmt", "yuv420p"])
ani.save(str(video_path), writer=writer, dpi=150)
plt.close(fig_anim)

print(f"動画保存: {video_path}")
print(f"\n出力ファイル:")
print(f"  静止画: {static_path}")
print(f"  動  画: {video_path}")
print(f"\n完了!")
