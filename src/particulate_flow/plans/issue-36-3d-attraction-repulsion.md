# Plan: issue-36 — 3D DEM 粒子間引力・斥力（Hamaker 型）

## Goal (verbatim from issue)

`DEM3D` (contact3d.py) の球-球間および球-円柱間に、2D `DEMSolver` と同じ Hamaker 型の
短距離引力・斥力を追加する。`LBMDEMSolver3D` と `runner3d.py` を通じて JSON config の
`attraction_strength` 等のパラメータが 3D 計算に有効になること。

## Hamaker 式（2D 準拠）

球-球間:
    r_eff = r_i * r_j / (r_i + r_j)   [接触距離 = r_i + r_j]
    h = max(surface_gap, max(attraction_min_gap, 1e-12))
    f = A* × r_eff / (6 × h²)         [スカラー大きさ]
    引力: 粒子 i を j に向けて引く（-n 方向）
    斥力: 粒子 i を j から遠ざける（+n 方向）

球-円柱間:
    r_eff = r_i * cyl_r / (r_i + cyl_r)
    同じ式。引力: 円柱表面に向かう方向（-normal_outward）

## Approach

### 1. `contact3d.DEM3D` に引数追加

```python
particle_attraction: bool = False,
particle_repulsion: bool = False,
attraction_strength: float = 1e-3,
repulsion_strength: float = 1e-3,
attraction_cutoff: float = 3.0,
repulsion_cutoff: float = 3.0,
attraction_min_gap: float = 0.05,
repulsion_min_gap: float = 0.05,
```

- `particle_attraction` と `particle_repulsion` 同時 True → `ValueError`（2D 準拠）
- これらを `self.` に保存

### 2. `DEM3D._sphere_sphere_loads` に引力・斥力項追加

既存の接触処理の *後*、`surface_gap = dist - (r_i + r_j)` を計算し、
2D `_apply_surface_force` 相当のロジックをインライン追加（contact3d.py 内に閉じる）。

surface_gap > cutoff の場合は接触がなくても引力が効く（接触前の近接域でも有効）ため、
ループは contact 判定の if 文の外で評価する。

### 3. `DEM3D._cylinder_loads` に引力・斥力項追加

球-円柱の `min_dist = cr + radii[i]` を使い、2D `_apply_cylinder_surface_force` 相当
のロジックをインライン追加。接触の有無によらず `surface_gap ≤ cutoff` で発動。

### 4. `lbm3d.LBMDEMSolver3D` への配線

`__init__` に同名の 8 引数を追加し、`DEM3D(...)` 呼び出しにそのまま渡す。
また `LBMDEMSolver3D` の docstring に引数の説明を追記する。
`particle_attraction` と `particle_repulsion` の排他チェックをここでも行う（lbm_dem.py 準拠）。

### 5. `runner3d.build_3d_solver` への配線

`getattr(args, "particle_attraction", False)` 等 8 個の `getattr` を追加。

## Alternatives Rejected

- **メソッド化**: 2D では `_apply_surface_force` / `_apply_cylinder_surface_force` として
  切り出しているが、3D は pure-Python ループ内のインラインが最小変更で済み、追加抽象化
  は過剰。
- **Numba カーネル拡張**: kernels.py の Numba カーネルに追加する案もあるが、今回は
  Pure-Python パスだけで十分（引力は O(N²) だが粒子数が少ない 3D ユースケース向け）。

## Out of Scope

- 粒子-壁面（y/z 境界面）引力（issue で明示的に除外）

## Test Plan

| Acceptance Scenario | テスト方法 |
|---|---|
| particle_attraction + strength 有効で引力が発生 | `DEM3D.compute_loads()` の力が attraction 方向を向く数値テスト |
| particle_repulsion 有効で斥力が発生 | 同上、逆向き |
| attraction と repulsion 同時 True → ValueError | `pytest.raises(ValueError)` |
| デフォルト (both False) で従来と同じ | 既存テストが通り続ける |
| config 経由で override できる | `runner3d.build_3d_solver` で `attraction_strength` が `dem.attraction_strength` に反映 |
| 球-円柱引力 | 円柱ありの `DEM3D` で引力方向の確認 |
