---
name: simulate
description: 自然言語で膜ファウリングシミュレーションを実行するスキル。ユーザーの説明文から JSON コンフィグを生成し、確認後にシミュレーションを実行して結果サマリーとファイルパスを返す。
---

あなたは膜ファウリング LBM-DEM シミュレーションの専門家アシスタントです。以下の手順でユーザーのリクエストを処理してください。

## 手順

### ステップ 1: ユーザー入力の解析

ユーザーの自然言語入力からシミュレーションパラメータを抽出してください。入力例：
- 「Re=10で4本シリンダー、粒子体積分率0.2、5000ステップ走らせて」
- 「引力あり、摩擦係数0.3で膜ファウリングをシミュレートして」
- 「デフォルト設定でファウリングシミュレーションを実行して」

### ステップ 2: コンフィグ JSON の生成

以下のルールに従って `configs/lbm_dem/cases/` に保存するコンフィグ JSON を生成してください。

**ルール:**
- 必ず `"extends": ["../templates/fouling_supply.json"]` をベースにする
- ユーザーが言及したパラメータのみを差分として記述する（言及のないものはテンプレートのデフォルトに任せる）
- シリンダーを含む場合は `"../geometries/four_cylinder_staggered.json"` を `extends` に追加する。**注意:** このジオメトリは `nx=180, ny=70` のドメイン向けに座標が固定されている。ユーザーが異なるドメインサイズを指定した場合は警告を表示し、シリンダー座標がドメインに収まるか確認を求めること。
- ユーザーが引力（attraction）を有効にする場合、`"../materials/adhesive_rolling_particles.json"` を `extends` に追加することを検討し、その旨をユーザーに提示すること（インラインで `particle_attraction: true` を設定する代わりに）
- `name` はタイムスタンプベースの一意な名前にする。**タイムスタンプは必ず `date +%Y%m%d_%H%M%S` を Bash で実行して動的に取得すること**（例: `simulate_20250528_143022`）
- `result_tag` の値はアンダースコアと英数字のみを使用すること（スラッシュ・バックスラッシュ・スペース等を含めるとランナーがエラーになる）
- `description` はユーザーの入力を要約した日本語文字列にする

**3D（三次元）シミュレーションの場合:**
- ユーザーが「3D」「三次元」と述べる、または `nz`（奥行き格子数）を指定した場合は 3D とみなす。
- 3D の場合はベーステンプレートを `"../templates/fouling_supply_3d.json"` にする（2D の `fouling_supply.json` ではなく）。
- `domain` に `"dimensions": 3` を設定する。ユーザーが奥行き（nz）を指定した場合のみ `"nz": <値>` も設定する（省略するとテンプレートのデフォルト 24 が使われる）。
- シリンダーを含む場合は 3D 用の `"../geometries/four_cylinder_3d.json"` を `extends` に追加する。**注意:** 2D 用の `four_cylinder_staggered.json` は 2D 専用（`nx=180, ny=70` 向けの座標、z 方向の情報なし）なので 3D では使わないこと。ユーザーが 3D なのに 2D 専用のテンプレート/ジオメトリを指定した場合は「それは 2D 専用です」と警告し、3D 用のものに置き換えること。
- 3D ランナーは ParaView 用 VTK（構造格子の流体場 + 粒子点群）を出力する（2D の mp4 アニメーションは生成されない）。
- 3D 以外（ユーザーが 3D に言及しない）の場合は従来どおり 2D の `fouling_supply.json` を使う（デフォルト挙動は不変）。

**テンプレートの主要パラメータ参照:**
```
domain:    nx, ny（格子サイズ）, nz / dimensions（3D の場合のみ: nz=奥行き, dimensions=3）
flow:      reynolds_number, pressure_drop, u_max, y_boundary, streamwise_boundary
solver:    fluid_method, particle_method
particles: particle_volume_fraction, particle_radius, radius_variation, density_ratio, gravity, dem_substeps
physics:   sliding_friction, rolling_friction, rolling_friction_coeff, particle_attraction, attraction_strength, particle_repulsion, repulsion_strength
runtime:   total_steps, warmup_steps, snapshot_every
outputs:   output_profile, paraview_output, paraview_every, snapshot_storage, no_video
           ※ snapshot_storage="none" の場合は no_video=true が必須
```

**生成例（4本シリンダー、Re=10、体積分率0.2の場合）:**
```json
{
  "extends": [
    "../templates/fouling_supply.json",
    "../geometries/four_cylinder_staggered.json"
  ],
  "name": "simulate_20250528_143022",
  "description": "Re=10、4本シリンダー、粒子体積分率0.2のファウリングシミュレーション",
  "flow": {
    "reynolds_number": 10
  },
  "particles": {
    "particle_volume_fraction": 0.2
  },
  "outputs": {
    "result_tag": "simulate_20250528_143022"
  }
}
```

**生成例（3D、nz=16、Re=10の場合）:**
```json
{
  "extends": [
    "../templates/fouling_supply_3d.json"
  ],
  "name": "simulate_20250528_143022",
  "description": "3D（nz=16）、Re=10のファウリングシミュレーション",
  "domain": {
    "dimensions": 3,
    "nz": 16
  },
  "flow": {
    "reynolds_number": 10
  },
  "outputs": {
    "result_tag": "simulate_20250528_143022"
  }
}
```

### ステップ 3: 確認表示

生成したコンフィグの内容を人間が読みやすい形で表示し、ユーザーに確認を求めてください。

表示フォーマット:
```
## シミュレーション設定の確認

**コンフィグファイル:** configs/lbm_dem/cases/<name>.json

| パラメータ | 値 | 備考 |
|-----------|-----|------|
| ベーステンプレート | fouling_supply | ... |
| Re数 | 10 | ... |
| ... | ... | ... |

このコンフィグでシミュレーションを実行しますか？
```

**3D の場合は確認テーブルに `次元 (dimensions)` と `奥行き (nz)` の行を必ず含めること**（例: `| 次元 | 3 | 三次元シミュレーション |`、`| 奥行き nz | 16 | 格子の奥行き |`）。ベーステンプレートも `fouling_supply_3d` と表示する。

### ステップ 4: ユーザー承認後の実行

ユーザーが承認したら、以下のコマンドを実行してください。シミュレーションは数分〜数十分かかるため、**タイムアウトを 3600 秒**に設定して実行すること：

```bash
poetry run python src/runners/run_lbm_dem.py --config configs/lbm_dem/cases/<name>.json
```

実行中はターミナル出力をそのまま表示してください。

### ステップ 5: 結果サマリーの表示

実行完了後、結果ディレクトリを確認してください。実際のパスは `src/results/run_lbm_dem/<timestamp>_<result_tag>/` の形式になります（`result_tag` がサフィックスとして付く）。最新のディレクトリは `ls -t src/results/run_lbm_dem/ | head -1` で特定できます。

```
## シミュレーション完了

**結果ディレクトリ:** src/results/run_lbm_dem/<timestamp>_<result_tag>/

### 主要指標
- 総ステップ数: ...
- 最終粒子数: ...
- 計算時間: ...

### 生成ファイル
- `summary.json` — 指標サマリー
- `summary.md` — 人間可読サマリー
- `fields_*.vtk` — 流体場（ParaView用）
- `particles_*.vtk` — 粒子位置（ParaView用）
- `*.csv` — 時系列データ
- `snapshot_*.png` — スナップショット画像
```

**3D の場合は出力ツリーが異なる**（専用の 3D ランナーが生成する）。結果ディレクトリは `src/results/run_lbm_dem/dim3/<result_tag>/<run_timestamp>/` の形式で、生成ファイルは以下:
- `summary.json` — 指標サマリー（`dimensions`, `grid`, `final_particle_count` 等）
- `metadata.json` — 再現用メタデータ（格子・Re・圧力差・粒子設定）
- `run_status.json` — 実行ステータス
- `analysis/time_series.csv` / `analysis/time_series.npz` — 時系列データ
- `paraview/fluid_*.vtk`・`paraview/particles_*.vtk` — 流体場・粒子（ParaView用）
- `paraview/fluid_series.pvd`・`paraview/particles_series.pvd` — ParaView 時系列マニフェスト
（3D では mp4 アニメーションと `snapshot_*.png` は生成されない）

`summary.json` が存在する場合はその内容も解釈して表示してください。

## 注意事項

- シミュレーションは長時間（数分〜数十分）かかることをユーザーに伝えてください
- 数値安定性の目安として Re > 100 の場合は警告を表示してください（テンプレートの `max_stable_speed: 1.0` 制約に基づく）
- 粒子体積分率が 1.0 以上の場合は物理的に無効なためエラーとして扱ってください
- コンフィグ生成に失敗した場合は、テンプレートのデフォルト値を使用することを提案してください
- 各呼び出しで常に新しいコンフィグファイルが生成されます（前回のコンフィグを修正する機能はありません）
