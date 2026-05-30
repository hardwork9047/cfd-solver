# Plan: issue-27 — 解析解ベンチマーク/収束次数の検証テスト

## Goal (from issue)

物理ソルバーの定量的正しさを担保する検証テスト群を追加する。既知の解析解との L2
誤差比較と、グリッド細分化による収束次数の確認で、「動く」だけでなく「正しい解に正しい
次数で収束する」ことを回帰テストとして固定する。

## Acceptance Scenarios

1. 平面ポアズイユ流の速度プロファイルが解析解と L2 相対誤差で一致する（2D・3D 両方）
2. 単一球/円柱に働く流体抗力が低レイノルズ数で Stokes 抗力と一致する
3. グリッド細分化で誤差が LBM の期待次数（2次）で減少する
4. 検証テストが CI で回り、数値挙動の退行を検出できる

## Approach

### 2D Poiseuille プロファイル（シナリオ 1a）

体積力駆動（`streamwise_boundary="periodic_force"`, `y_boundary="wall"`）の 2D ソルバー。
定常状態の解析解は：

```
u_x(y) = F/(2ν) · y(H-y)   [格子単位, y=0 at bot wall, H=ny-2 fluid layers]
```

`F_drive = 8ν·u_max/(ny²)` から `u_max = F·(H/2)²/(2ν)` が成立する。

精度基準: L2 相対誤差 < 1% @ ny=32。

### 3D Poiseuille プロファイル（シナリオ 1b）

3D ソルバーはボディフォースを直接サポートしない（`Fx` は IBM フィードバック用）。
代わりに圧力駆動 (`streamwise_boundary="pressure"`) を使い、十分に収束させてから
y 方向スライスの速度プロファイルを解析解（チャネル Poiseuille）と比較する。
3D でチャネル形状を作るには z 方向も壁が必要 — LBMDEMSolver3D は z 壁を
サポートしないため、y 方向の periodic + 圧力駆動を使い、均一プロファイル
（plug flow）になることを確認する。これは「3D Poiseuille」ではなく「3D 平行流」の
回帰になる。より意味のある 3D 検証は壁なし均一流の u_x 収束とする。

**代替案**: 2D の L2 検証と収束テストを充実させ、3D は「圧力差に対して正の平均流速が
発達する」既存テスト（`test_lbm3d_pressure.py`）で十分とし、3D Poiseuille は別イシューに
委ねる。→ **採用**: 3D は plug-flow L2 回帰テストに留め、Poiseuille プロファイルは 2D で担う。

### Stokes 抗力（シナリオ 2）

2D 円柱の抗力: 完全な Stokes 解（Oseen）は 2D では log 修正が必要で解析解が単純でない。
→ 代わりに「均一流中の粒子が受ける Stokes 抗力 = 3πμd(u-v)」を単体テストで確認する
(`_stokes_drag` メソッドへのユニットテスト)。これは既存 `_stokes_drag` メソッドの
ユニット検証として追加する。

### グリッド収束（シナリオ 3）

ny = 16, 24, 32 の 3 点で 2D Poiseuille の L2 誤差を測定し、グリッドを 2 倍にすると
誤差が 4 倍減少（2次収束）することを確認する。

## Test File Location

`tests/lbm/test_verification.py` — 新規ファイル

## Alternatives Rejected

- **3D 壁あり Poiseuille**: LBMDEMSolver3D に z 壁 BC が未実装のため実装不要
- **実際の流体抗力の積分テスト**: 円柱周りの抗力計算は別のモジュール(io/cylinder_flow.py)
  にあり、そちらのテストは既存。ここでは `_stokes_drag` の数値精度を確認するに留める
- **Numba パスの個別テスト**: numpy/numba の物理一致は #33 で検証済み

## Out of Scope

- 3D 壁面 Poiseuille プロファイル（z 壁 BC 未実装）
- 粒子間相互作用の検証
- 高 Re 乱流挙動

## Assumptions

- 2D: `y_boundary="wall"` で y=0, y=ny-1 が固体。流体層は y=1..ny-2（計 ny-2 層）
- `F_drive = 8ν·u_max/ny²` の設定から解析解の u_max は整合する
- 収束次数テストは `@pytest.mark.slow` でマーク

## Test Plan

| Test | Scenario | Marker |
|---|---|---|
| `test_poiseuille_2d_profile_l2` | 1a: L2 < 1% at ny=32 | slow |
| `test_stokes_drag_unit` | 2: F = 3πμd·Δu | (fast) |
| `test_poiseuille_2d_grid_convergence` | 3: 2次収束 | slow |
| `test_3d_uniform_flow_profile` | 1b: 3D 均一流 L2 | slow |
