# Plan: Issue-1 — particulate_flow サブパッケージ分割

## Goal

`particulate_flow/` を責務別サブパッケージに分割し、AI がコードを読む際のトークン消費を削減する。
冗長なインラインコメントを整理し、テストをサブディレクトリ構成に移行する。
全既存テスト・公開API・NumPy/Numba切り替えを維持する。

---

## 現状の問題

- `lbm_dem.py` が 2,363 行・64 メソッドのモノリシッククラス
- D2Q9定数、Numbaカーネル、LBM衝突、IBM結合、粒子管理、可視化が1ファイルに混在
- AI が任意の機能を読む際に全体コンテキストが必要になる

---

## ターゲット構成

```
src/particulate_flow/
  __init__.py              # 既存公開API を re-export（後方互換）
  solver.py                # LBMDEMSolver（オーケストレータ）
  fast_solver.py           # FastLBMDEM（変更なし）
  lbm/
    __init__.py
    constants.py           # C, W, OPPOSITE, Q, CS2, FLUID_METHODS 等
    operators.py           # _equilibrium, _guo_forcing（純関数）
    kernels.py             # Numba LBM カーネル（_lbm_step_numba）
    collision.py           # _collide_bgk_guo, _collide_trt_guo, _lbm_step
    boundary.py            # _apply_pressure_boundaries
    macroscopic.py         # _macroscopic, _fluid_max_speed, _control_drive_force
  dem/
    __init__.py
    solver.py              # DEMSolver（dem_solver.py から移動）
    packing.py             # DEMPackingSimulation（dem_packing.py から移動）
    contact.py             # DEM 接触・セルリスト・_dem_substep 関連メソッド群
    particle_manager.py    # 粒子の追加・削除・inlet管理メソッド群
    kernels.py             # Numba DEM カーネル（_dem_pair_loads_numba）
  ibm/
    __init__.py
    coupling.py            # IBM 結合ロジック（_apply_immersed_boundary_forces 等）
    kernels.py             # Numba IBM カーネル（_ibm_markers_numba 等）
  geometry/
    __init__.py
    pore.py                # PoreGeometry, Cylinder（geometry.py から移動）
  io/
    __init__.py
    config.py              # SimulationConfig（simulation_config.py から移動）
    paths.py               # result_paths（result_paths.py から移動）
    verification.py        # fluid_verification（fluid_verification.py から移動）
    visualization.py       # plot_fields, plot_particles, cylinder_flow
```

---

## アプローチ

### フェーズA: 既存ファイルをサブパッケージへ移動

移動対象（実装変更なし、import パス変更のみ）:
| 現在 | 移動先 |
|------|--------|
| `dem_solver.py` | `dem/solver.py` |
| `dem_packing.py` | `dem/packing.py` |
| `geometry.py` | `geometry/pore.py` |
| `simulation_config.py` | `io/config.py` |
| `result_paths.py` | `io/paths.py` |
| `fluid_verification.py` | `io/verification.py` |
| `cylinder_flow.py` | `io/visualization.py` に統合 |

### フェーズB: `lbm_dem.py` の分割

モジュールレベルのコードを抽出:
- D2Q9定数 → `lbm/constants.py`
- `_equilibrium`, `_guo_forcing` → `lbm/operators.py`
- Numba LBM カーネル → `lbm/kernels.py`
- Numba DEM カーネル → `dem/kernels.py`
- Numba IBM カーネル → `ibm/kernels.py`
- `plot_fields`, `plot_particles` 等 → `io/visualization.py`

`LBMDEMSolver` のメソッドをミックスインクラスとして抽出:
- LBM衝突・境界・マクロ量 → `lbm/` 内のミックスインまたはヘルパーモジュール
- DEM接触・粒子管理 → `dem/contact.py`, `dem/particle_manager.py` 内のスタンドアロン関数
- IBM結合 → `ibm/coupling.py` 内のスタンドアロン関数

`LBMDEMSolver` 自体は `solver.py` に残し、各サブパッケージの関数を `import` して呼び出すオーケストレータとして機能させる。

### フェーズC: コメント整理

対象: 意図が型ヒント・変数名から自明なコメントを削除または短縮。
除外: docstring、物理的な意味の説明、参考文献。

### フェーズD: テスト移行

```
tests/
  conftest.py
  lbm/
    __init__.py
    test_collision.py      # 衝突・境界テスト（test_lbm_dem.py から分割）
    test_macroscopic.py
  dem/
    __init__.py
    test_solver.py         # test_lbm_dem.py（DEM部）+ test_dem_packing.py
    test_contact.py
  ibm/
    __init__.py
    test_coupling.py
  geometry/
    __init__.py
    test_pore.py           # test_geometry.py から移動
  io/
    __init__.py
    test_config.py         # test_simulation_config.py から移動
  test_integration.py      # エンドツーエンド統合テスト（既存の slow テスト）
```

---

## 代替案と不採用理由

**LBMDEMSolver を複数クラスに分解する**
→ フェーズ2（runner整理）で依存するインターフェースが変わるため、今は行わない。
オーケストレータとして維持し、フェーズ2でコンポーネント化する。

**フラット構成のままファイルを増やす**
→ サブパッケージなしではディレクトリが散らかり、AI が「この機能はどこ？」と推測しなければならない。

---

## スコープ外

- runner の変更（フェーズ2）
- configs の変更（フェーズ3）
- `LBMDEMSolver` の公開 API 変更（後方互換を維持）
- パフォーマンス最適化

---

## テスト戦略

各 Acceptance Scenario に対応:
1. `from particulate_flow.lbm import ...` → 各サブパッケージの `__init__.py` に import テスト
2. 後方互換 → `from particulate_flow import LBMDEMSolver` が動作することを既存テストで確認
3. NumPy/Numba 切り替え → `--fluid-accelerator numba` / `numpy` のスモークテスト
4. コメント整理 → 主観的なので自動テスト対象外（コードレビューで確認）
5. テスト移行 → `pytest tests/` で全テストがパスすること

---

## 実装順序

1. サブパッケージディレクトリと `__init__.py` を作成
2. 既存ファイルをサブパッケージへ移動（フェーズA）
3. `lbm_dem.py` からモジュールレベルコードを抽出（フェーズB前半）
4. `LBMDEMSolver` のメソッドをヘルパー関数として各サブパッケージに抽出（フェーズB後半）
5. `__init__.py` で後方互換 re-export を設定
6. テストを移行（フェーズD）
7. コメント整理（フェーズC）
