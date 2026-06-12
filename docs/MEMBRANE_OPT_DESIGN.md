# membrane_opt プログラム設計書

[`TOPOLOGY_OPTIMIZATION_PLAN.md`](TOPOLOGY_OPTIMIZATION_PLAN.md)(研究計画)を実装するためのプログラム設計。
src 配下の全モジュールを精査した結果に基づく。

---

## 1. 設計を支えるコードベースの事実

調査で確認した、設計判断の根拠となる4点:

1. **形状はソルバー改変ゼロで注入できる**
   config の `geometry.cylinders` ブロックは `cylinder_spec` に正規化され
   (`src/particulate_flow/io/config.py` の `_normalise_geometry_cylinders`)、
   `extends` 継承も解決される。形状生成器は「ケース JSON を吐くだけ」でよい。

2. **乱数シードが固定されている(唯一のコア側修正点)**
   - 2D: `src/particulate_flow/builder.py` の `build_lbm_dem_solver` 内で `seed=42` ハードコード
   - 3D: `LBMDEMSolver3D(seed=42)` デフォルトのまま `runner3d.build_3d_solver` が渡していない

   Phase 0 のノイズ測定(同一形状・シード違い反復)には `--seed` の公開が必須。

3. **幾何記述子解析の原型が既に存在する**
   `src/tools/analyze_lbm_dem_design_sweeps.py` の `_geometry_features`
   (円柱数・平均半径・固体面積率・最小ギャップ)と Spearman 相関ランキング。
   層2(設計則抽出)はこの移植・拡張で済む。

4. **実行系は2系統を使い分けられる**
   - 安価な 2D DoE バッチ: `src/runners/run_lbm_dem_sweep.py --jobs N`(単一マシン並列)
   - 高価な 3D 検証: `src/runners/sweep_cli.py` マニフェスト(複数マシン git 分担)

   最適化ループはどちらにも「スイープ定義/マニフェストを書き出す」だけで接続できる。

その他の確認事項:
- `removed_particles` は「出口を完全に通過した粒子の削除数」(`runner3d.py` の
  `_remove_outflow_particles`)= **透過粒子数**。阻止率の分母分子はここから取る。
- 円柱が圧力境界(x=0, x=nx-1)に達するとソルバーがエラーを出す(`lbm3d.py` の
  `_build_cylinder_solid` 検証)→ 設計空間の制約として先取りする。

---

## 2. 全体構造 — 新パッケージ `src/membrane_opt/` + ランナー 1 本

最適化はシミュレータの「消費者」なので `particulate_flow` には入れず独立パッケージとする
(pyproject の `packages` に `{include = "membrane_opt", from = "src"}` を追加)。
Web UI が sweep CLI を薄くラップしたのと同じ思想で、**BO ループは既存スイープ基盤を薄くラップ**する。

```
src/membrane_opt/
  design.py       設計空間の定義・エンコード/デコード・制約検証
  generate.py     設計ベクトル → ケース JSON・マニフェスト生成
  metrics.py      結果ディレクトリ → 目的関数 (f1, f2, f3)
  descriptors.py  幾何記述子(層2)
  surrogate.py    GP + 多目的 BO(層1、botorch)
  ledger.py       スタディ台帳(評価履歴の追記型 CSV)
src/runners/run_membrane_opt.py   オーケストレータ CLI
configs/membrane_opt/             スタディ定義・基底テンプレート
optimization/<study>/             台帳・提案・パレート出力(git 追跡)
```

## 3. データフロー

```
study.json(設計空間・目的・予算・忠実度)
   │
   ▼ suggest                          ▼ 実行(既存基盤・改変なし)
BO/LHS が設計バッチを提案 ──→ ケースJSON + マニフェスト ──→ 2D: run_lbm_dem_sweep.py --jobs N
   ▲                                                        3D: ./bin/sweep run(複数マシン)
   │ harvest                                                      │
ledger.csv に (設計ベクトル, seed, f1, f2, f3) を追記 ←── time_series.csv / summary.json
   │
   ▼ pareto
パレートフロント CSV + 散布図 + 代表形状の geometry JSON
```

---

## 4. モジュール別の責務

### `design.py` — 設計空間

- `DesignSpace`: 変数の bounds と制約を持つデータクラス
  - Stage A: 円柱 N 本の (x, y, r)(10〜30 次元)
  - Stage B: 構造化パラメータ(千鳥度・列数・半径分布・のど幅プロファイル)
- `decode(x) -> list[Cylinder]`: 既存 `src/particulate_flow/geometry/pore.py` の
  `Cylinder` / `PoreGeometry` を再利用
- `validate(x) -> list[str]`: 制約違反の列挙
  - 円柱非重複
  - 最小のど幅 ≥ 粒子径 × 係数(粒子が物理的に通れる下限)
  - 空隙率レンジ
  - 入出口クリアランス(圧力境界エラーの先取り)

### `generate.py` — 設計 → ケース config

- `write_case(design, seed) -> Path`: 基底テンプレートへの `extends` +
  `geometry.cylinders` ブロック + `seed` だけの薄い JSON を生成
- ケース名 = `opt_<study>_<設計ハッシュ>_s<seed>`(決定論的 → 再実行・重複検出が容易)
- 出力先を 2 系統サポート:
  - 2D: `run_lbm_dem_sweep.py` 用スイープ JSON(`cases` リスト)
  - 3D: `sweep_cli` マニフェスト(`configs/lbm_dem/sweeps/manifests/`)

### `metrics.py` — 結果 → 目的関数

`evaluate(case_dir) -> {f1, f2, f3}`(+ 反復集約)

| 目的 | 算出 | データ源 |
|---|---|---|
| f1 初期透過性 | ウォームアップ後初期の passed_rate = `diff(removed_particles)/diff(step)` | `analysis/time_series.csv` |
| f2 阻止率 | `1 − removed_particles / generated_particles` | `summary.json` |
| f3 透過性維持率 | 後期 passed_rate / 初期 passed_rate | `analysis/time_series.csv` |

- 同一設計の複数 seed を平均・分散付きで集約(noisy GP の入力にする)

### `descriptors.py` — 幾何記述子(層2)

- `analyze_lbm_dem_design_sweeps.py` の `_geometry_features` を移植・拡張:
  空隙率、最小のど幅/粒径比、千鳥度、迂曲度プロキシ、固体周長
- 用途は 2 つ: 設計則回帰(scikit-learn)の特徴量、BO の補助入力

### `surrogate.py` — サロゲート + 提案(層1)

- botorch の **qNEHVI**(ノイズ対応・多目的・バッチ提案)
- 共通インターフェース `Sampler.propose(q) -> list[DesignVector]` に
  LHS サンプラー(Phase 1 DoE 用)と BO サンプラーを両方載せる
  → DoE → BO の切り替えは study.json の 1 キー

### `ledger.py` — スタディ台帳

- `optimization/<study>/ledger.csv` への追記型台帳
  (設計ベクトル、seed、目的値、ケース名、状態)
- `sweep_results/` と同じ「git 同期される軽量状態」哲学
  → 複数マシンでの 3D 評価でも台帳をマージできる

### `run_membrane_opt.py` — オーケストレータ CLI

| サブコマンド | 動作 |
|---|---|
| `init` | スタディ作成(設計空間・目的・予算・忠実度を study.json に) |
| `suggest` | サンプラーが設計バッチを提案 → ケース JSON + スイープ定義を書き出し |
| `harvest` | 完了ケースの結果を f1/f2/f3 に変換して台帳へ追記 |
| `status` | 台帳の進捗表示 |
| `pareto` | パレートフロント CSV + 散布図 + 代表形状の geometry JSON を出力 |
| `auto` | 2D 限定: suggest → 実行 → harvest を自動ループ |

3D は `suggest` で提案だけ出し、実行は人間が `./bin/sweep run` で分担する半自動運用。

---

## 5. コア側の変更(最小限・1 点のみ)

`--seed`(default 42)を公開する:

1. `src/runners/run_lbm_dem.py` に `--seed` 引数を追加
2. `src/particulate_flow/builder.py` のハードコード `seed=42` を `getattr(args, "seed", 42)` に
3. `src/particulate_flow/runner3d.py` の `build_3d_solver` で `LBMDEMSolver3D(seed=...)` に渡す
4. `metadata.json` に seed を記録

約 10 行、後方互換(デフォルト 42 で既存結果と同一挙動)。

---

## 6. 依存関係

web グループと同じパターンで optional group を追加:

```toml
[tool.poetry.group.opt]
optional = true

[tool.poetry.group.opt.dependencies]
botorch = ">=0.11"      # torch を引き込む(~2GB)。BO を使うマシンのみ install --with opt
scikit-learn = ">=1.4"  # 層2の記述子回帰
```

**注:** Phase 0–1(seed 反復・LHS・記述子回帰)は botorch 不要。
torch の重さが気になる場合は導入を Phase 2 まで遅延できる
(scikit-learn のみの軽量グループから始める選択肢もある)。

---

## 7. フェーズ対応表

| フェーズ | 実装するもの | 再利用する既存資産 |
|---|---|---|
| Phase 0 | `--seed` 公開、`metrics.py`、`generate.py`(反復ケース生成) | sweep 基盤で実行 |
| Phase 1 | `design.py`、`descriptors.py`、LHS サンプラー | `run_lbm_dem_sweep.py --jobs`、`_geometry_features` |
| Phase 2 | `surrogate.py`(qNEHVI)、`auto` ループ、`pareto` | `sweep_cli` マニフェスト(3D 検証) |
| Phase 3 | ボクセルマスク形状 + CNN サロゲート | `lbm3d.py` `_build_cylinder_solid` の拡張 |

---

## 8. テスト方針

実 sim を起動しないテストで固める(`tests/test_sweep_web.py` と同じ流儀):

- 設計ベクトル ⇔ 円柱リストのラウンドトリップ
- 制約検証(重複・のど幅・クリアランス違反の検出)
- 生成ケースが `SimulationConfig.from_json` で解決できること
- 合成 time_series.csv からの f1/f2/f3 算出
- 台帳の追記・重複検出・マージ
- BO は解析的トイ関数(ZDT 系)での煙テスト(botorch 未導入環境では `importorskip`)

---

## 9. 着手順

実装に入る場合の最初の一歩は **Phase 0**:

1. `--seed` 公開(コア側 ~10 行)
2. `metrics.py` + テスト
3. `generate.py` で baseline 形状 × seed 5 反復のマニフェスト生成
4. sweep 基盤で実行 → ノイズ分散を測定 → GP のノイズ項と必要反復数を決定

---

## 関連ファイル

| パス | 内容 |
|---|---|
| `docs/TOPOLOGY_OPTIMIZATION_PLAN.md` | 研究計画(本書はその実装設計) |
| `docs/RESEARCH_ROADMAP.md` | メカニズム解明ロードマップ |
| `src/particulate_flow/io/config.py` | config 正規化・extends 解決(生成器の出力先) |
| `src/particulate_flow/geometry/pore.py` | `Cylinder` / `PoreGeometry`(design.py が再利用) |
| `src/particulate_flow/builder.py` | seed ハードコード箇所(Phase 0 修正対象) |
| `src/particulate_flow/runner3d.py` | 3D 実行パス・seed 未伝搬(Phase 0 修正対象) |
| `src/tools/analyze_lbm_dem_design_sweeps.py` | `_geometry_features`(descriptors.py の母体) |
| `src/runners/run_lbm_dem_sweep.py` | 2D DoE バッチ実行系 |
| `src/runners/sweep_cli.py` | 3D 複数マシン実行系 |
