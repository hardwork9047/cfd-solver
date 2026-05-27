# Plan: Issue-3 — src/bin/ 廃止・src/tools/ 新設・configs/ 再整理

## Goal

`src/bin/` を廃止し、config sweep runner を `src/runners/` へ、分析・ベンチマーク系スクリプトを
`src/tools/` へ移動する。`configs/` を `cases/`, `geometries/`, `materials/`, `templates/`,
`sweeps/` の共通リソース構造に再整理する。

---

## 現状分析

### src/bin/ (11 Python + 3 shell スクリプト)

| ファイル | 分類 | 移動先 |
|---|---|---|
| `run_lbm_dem_config_sweep.py` | sweep runner | `src/runners/run_lbm_dem_sweep.py` |
| `analyze_lbm_dem_design_sweeps.py` | 分析 | `src/tools/analyze_lbm_dem_design_sweeps.py` |
| `benchmark_ibm_marker_spacing.py` | ベンチ | `src/tools/benchmark_ibm_marker_spacing.py` |
| `benchmark_lbm_accelerators.py` | ベンチ | `src/tools/benchmark_lbm_accelerators.py` |
| `benchmark_lbm_dem_coupling_backends.py` | ベンチ | `src/tools/benchmark_lbm_dem_coupling_backends.py` |
| `benchmark_numba_components.py` | ベンチ | `src/tools/benchmark_numba_components.py` |
| `plot_lbm_dem_fouling_metrics.py` | 分析 | `src/tools/plot_lbm_dem_fouling_metrics.py` |
| `summarize_lbm_dem_results.py` | 分析 | `src/tools/summarize_lbm_dem_results.py` |
| `verify_lbm_dem_fluid.py` | 検証 | `src/tools/verify_lbm_dem_fluid.py` |
| `run_dem_packing_200.sh` | shell runner | 削除（configs+runners で代替） |
| `run_lbm_dem_membrane_pressure_periodic.sh` | shell runner | 削除（smoke configに代替あり） |
| `run_lbm_dem_solver_method_examples.sh` | shell runner | 削除 |
| `run_particulate_flow_cylinder_flow.sh` | shell runner | 削除 |

### configs/ 現状

```
configs/
  lbm_dem/
    cases/       cylinder_flow_single.json, fouling_four_cylinder_supply.json, dem_settling_200.json
    geometries/  four_cylinder_staggered.json
    materials/   adhesive_rolling_particles.json
    templates/   cylinder_flow.json, dem_settling_pack.json, fouling_supply.json
    sweeps/      fouling_screening_example.json
    membrane_pressure_periodic_smoke.json   ← 分類されていないケース
  dem_packing/
    cases/       packing_200_particles.json
```

問題点:
- `configs/lbm_dem/membrane_pressure_periodic_smoke.json` が `cases/` に分類されていない
- `dem_packing/` 以下に `templates/`, `geometries/` 等がない（現時点では1ケースのみ、追加不要）
- 全体の `configs/cases/`, `configs/templates/` への統合は不採用（タイプ別サブディレクトリを維持）

### extends 参照パスの影響

変更前 → 変更後:
- `cases/fouling_four_cylinder_supply.json` の `"../templates/fouling_supply.json"` → 変更なし（lbm_dem内の相対パスは維持）
- `membrane_pressure_periodic_smoke.json` が `cases/` に移動 → 参照側なし（この設定は他から extends されていない）

---

## アプローチ

### 1. `src/tools/` 新設、Python スクリプト移動

`src/tools/__init__.py` なし（スタンドアロンスクリプト群）。
各スクリプトの `sys.path.insert(0, ...)` の `parents[N]` インデックスを調整する。

### 2. `src/runners/run_lbm_dem_sweep.py` 追加

`src/bin/run_lbm_dem_config_sweep.py` をそのまま移動。
`REPO_ROOT = Path(__file__).resolve().parents[2]` → `.parents[2]` 変わらず（src/runners/ から2上がれば repo root）。

### 3. `configs/lbm_dem/membrane_pressure_periodic_smoke.json` を `cases/` に移動

ファイル移動のみ。中身変更なし。

### 4. `src/bin/` 削除

shell スクリプト含め `src/bin/` ディレクトリごと削除。

### 5. CLAUDE.md 更新

`src/bin/` への参照を削除・更新。`src/tools/` の記述追加。
`configs/lbm_dem/membrane_pressure_periodic_smoke.json` のパスを更新。

---

## 代替案と不採用理由

**configs/ を `configs/cases/`, `configs/templates/` と全タイプ横断で統合する**
→ `lbm_dem/` と `dem_packing/` でテンプレートの形式が異なるため、タイプ別を維持する方が明確。
  現状でも `lbm_dem/cases/`, `lbm_dem/templates/` 等が既にある。

**shell スクリプトを runners に移動する**
→ shell スクリプトは特定の引数を固定して呼び出す例として機能していたが、JSON configs で
  同等の再現が可能なため削除が適切。

---

## スコープ外

- `configs/` 内の値の変更（内容は変えない）
- `dem_packing/` への `templates/`, `geometries/` 追加（ケースが1つのみ、不要）
- `src/tools/` スクリプトへの大きな機能追加

---

## テスト戦略

各 Acceptance Scenario に対応:

1. `src/bin/` が存在しない → テストでディレクトリ存在確認
2. `src/tools/` のスクリプトが import error なし → 各スクリプトの compile check
3. `configs/lbm_dem/cases/` に smoke config が含まれる → ファイル存在確認
4. `configs/geometries/` 等が独立して存在する → ディレクトリ存在確認（既存の lbm_dem/ 内）
5. `"extends"` 参照パスが正しい → `SimulationConfig.from_json()` で各 cases を読み込めるか確認
6. CLAUDE.md にパスが記述されている → テスト対象外（manual review）

---

## 実装順序

1. `src/tools/` 新設、Python スクリプト9本移動（sys.path 調整）
2. `src/runners/run_lbm_dem_sweep.py` 追加（sweep runner 移動）
3. `configs/lbm_dem/membrane_pressure_periodic_smoke.json` を `cases/` に移動
4. `src/bin/` 削除
5. CLAUDE.md 更新
6. テスト追加・確認
