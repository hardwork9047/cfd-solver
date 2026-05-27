# Plan: Issue-2 — runner エントリポイント化と組み合わせ対応

## Goal

`src/runners/` の runner をエントリポイントのみの薄いスクリプトにし、ソルバーの
組み立てロジックを `particulate_flow/` 内に移す。
JSONコンフィグに `solver`、`accelerator` セクションを追加し、モデル・数値手法・
geometry を設定で組み合わせられるようにする。

---

## 現状分析

| ファイル | 行数 | 主な内容 |
|---------|------|---------|
| `src/runners/run_lbm_dem.py` | 1,769行 | CLI引数定義・PoreGeometry構築・FastLBMDEM構築・出力/分析/アニメーション |
| `src/runners/run_dem_packing.py` | 270行 | CLI引数定義・DEMPackingSimulation構築・出力 |
| `src/particulate_flow/io/config.py` | 201行 | `SimulationConfig` + `extends` 継承 |

**現在の config `SECTION_KEYS`**: domain, flow, numerics, particles, physics, runtime, outputs, stability

`numerics` セクションに `fluid_method`, `fluid_accelerator`, `compute_accelerator`,
`particle_method`, `particle_fluid_coupling` が混在している。

---

## アプローチ

### 1. `SimulationConfig` に `solver` / `accelerator` セクションを追加

**`solver` セクション** → 物理モデルの選択:
```json
"solver": {
  "fluid_method": "lbm-trt-guo",
  "particle_method": "dem-hertz",
  "particle_fluid_coupling": "immersed_boundary",
  "ibm_marker_spacing": 2.0
}
```

**`accelerator` セクション** → バックエンドの選択:
```json
"accelerator": {
  "fluid": "auto",
  "compute": "auto"
}
```

`accelerator` セクションのキー `fluid` → `fluid_accelerator`、`compute` → `compute_accelerator`
に `_normalise_key()` でマッピングする。

既存の `numerics` セクションも後方互換で維持する（`fluid_method` 等は `numerics` にあっても動く）。

### 2. `particulate_flow/builder.py` を新規作成

`SolverConfig` dataclass と `build_lbm_dem_solver()` 関数を提供する:
- `args` (argparse.Namespace) から `PoreGeometry` を構築
- `FastLBMDEM(...)` を構築して返す
- runner はこれを呼ぶだけ

同様に `build_dem_packing_solver()` も追加。

### 3. runner をエントリポイント化

`run_lbm_dem.py`:
- `sim = FastLBMDEM(nx=NX, ...)` の 50行 → `sim = build_lbm_dem_solver(args)` 1行
- `GEOMETRY = PoreGeometry.from_cylinders(CYLINDERS)` → builder 内に移動
- 出力/分析/アニメーション (~1700行) はそのまま維持

`run_dem_packing.py`:
- `DEMPackingSimulation(...)` の構築 → `build_dem_packing_solver(args)` に移動

### 4. 既存テンプレートに `solver` / `accelerator` セクションを追加

`configs/lbm_dem/templates/fouling_supply.json` 等に `solver` / `accelerator` セクションを追加。
`numerics` から solver 系キーを `solver` セクションへ移動。
後方互換のため `numerics` に残ったキーは引き続き動作する。

---

## 代替案と不採用理由

**runner の出力・分析コードも `particulate_flow/io/` に移動する**
→ runner が 1,700行超のアニメーション・VTK・CSV出力ロジックを持つが、これは
  出力フォーマット固有のコードであり、フェーズ2のスコープを超える。フェーズ3で対応。

**`solver` セクションを `numerics` にマージする**
→ 物理モデル選択（何を計算するか）とバックエンド選択（どう計算するか）の関心が混在するため分離する。

---

## スコープ外

- `src/bin/` の整理（フェーズ3）
- `configs/` 全体の再整理（フェーズ3）
- runner の出力/分析コードの `particulate_flow/` への移動
- runner の行数を大幅削減すること（出力コードは維持）

---

## テスト戦略

各 Acceptance Scenario に対応:

1. `run_lbm_dem.py` の組み立てロジックなし → `build_lbm_dem_solver()` 関数のユニットテスト
2. `solver.fluid_method` でBGK/TRT切り替え → `SimulationConfig.from_mapping({"solver": {"fluid_method": "lbm-bgk-guo"}})` テスト
3. `accelerator.fluid` でバックエンド切り替え → `SimulationConfig.from_mapping({"accelerator": {"fluid": "numpy"}})` テスト
4. `geometry` セクション → 既存テストを活用 + builder テスト
5. `"extends"` 継承 → 既存の `test_simulation_config` テストで確認
6. サブパッケージ利用 → `build_lbm_dem_solver` が `particulate_flow.lbm`, `particulate_flow.dem` 等を import するテスト

---

## 実装順序

1. `SimulationConfig` に `solver` / `accelerator` セクション追加
2. `particulate_flow/builder.py` 作成
3. `run_lbm_dem.py` で `build_lbm_dem_solver()` を使用
4. `run_dem_packing.py` で `build_dem_packing_solver()` を使用
5. テンプレートJSONに新セクション追加
