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
- シリンダーを含む場合は `"../geometries/four_cylinder_staggered.json"` を `extends` に追加する
- `name` はタイムスタンプベースの一意な名前にする（例: `simulate_20240528_143022`）
- `description` はユーザーの入力を要約した日本語文字列にする

**テンプレートの主要パラメータ参照:**
```
domain:    nx, ny（格子サイズ）
flow:      reynolds_number, pressure_drop, u_max, y_boundary, streamwise_boundary
solver:    fluid_method, particle_method
particles: particle_volume_fraction, particle_radius, radius_variation, density_ratio, gravity
physics:   sliding_friction, rolling_friction, rolling_friction_coeff, particle_attraction, attraction_strength, particle_repulsion, repulsion_strength
runtime:   total_steps, warmup_steps, snapshot_every
outputs:   output_profile, paraview_output, no_video
```

**生成例（4本シリンダー、Re=10、体積分率0.2の場合）:**
```json
{
  "extends": [
    "../templates/fouling_supply.json",
    "../geometries/four_cylinder_staggered.json"
  ],
  "name": "simulate_20240528_143022",
  "description": "Re=10、4本シリンダー、粒子体積分率0.2のファウリングシミュレーション",
  "flow": {
    "reynolds_number": 10
  },
  "particles": {
    "particle_volume_fraction": 0.2
  },
  "outputs": {
    "result_tag": "simulate_20240528_143022"
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

### ステップ 4: ユーザー承認後の実行

ユーザーが承認したら、以下のコマンドを実行してください：

```bash
poetry run python src/runners/run_lbm_dem.py --config configs/lbm_dem/cases/<name>.json
```

実行中はターミナル出力をそのまま表示してください。

### ステップ 5: 結果サマリーの表示

実行完了後、`src/results/run_lbm_dem/<timestamp>/` に生成された結果を確認し、以下の形式でサマリーを表示してください：

```
## シミュレーション完了

**結果ディレクトリ:** src/results/run_lbm_dem/<timestamp>/

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

`summary.json` が存在する場合はその内容も解釈して表示してください。

## 注意事項

- シミュレーションは長時間（数分〜数十分）かかることをユーザーに伝えてください
- 数値的に不安定になりそうな設定（Re数が高すぎる、体積分率が1.0以上など）は警告を表示してください
- コンフィグ生成に失敗した場合は、テンプレートのデフォルト値を使用することを提案してください
