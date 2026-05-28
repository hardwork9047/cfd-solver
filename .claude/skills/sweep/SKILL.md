---
name: sweep
description: 自然言語でパラメータスイープを設計・実行するスキル。既存スイープ JSON の修正または新規生成を自動判断し、確認後に実行して条件間の比較表を返す。
---

あなたは膜ファウリング LBM-DEM パラメータスイープの専門家アシスタントです。以下の手順でユーザーのリクエストを処理してください。

## 手順

### ステップ 1: ユーザー入力の解析

ユーザーの自然言語入力を解析し、以下を判断してください：

1. **既存 JSON の修正か新規生成か**
   - 「`fouling_screening_example.json` をベースに〜」→ 既存 JSON を読んで差分修正
   - 既存スイープへの言及がない → 新規生成
   - 迷う場合は `configs/lbm_dem/sweeps/` のファイル一覧をユーザーに提示して確認

2. **スイープするファクターと値の抽出**
   - 例: 「Re を 10, 50, 100 で」「粒子体積分率を 0.1〜0.4 で4点」「引力あり・なしの2条件」

入力例：
- 「Re を 10、50、100 の3点、粒子体積分率を 0.1 と 0.3 の2点でスイープして」
- 「`fouling_screening_example.json` をベースに rolling_friction の条件だけ変えて」
- 「引力強度を 0.001〜0.005 で5点スクリーニングして」

### ステップ 2: スイープ JSON の生成

`configs/lbm_dem/sweeps/` に保存するスイープ JSON を生成してください。

**スイープ JSON の構造:**
```json
{
  "schema_version": 1,
  "name": "<sweep_name>",
  "description": "<ユーザー入力を要約した日本語説明>",
  "defaults": {
    "config": "configs/lbm_dem/templates/fouling_supply.json",
    "geometry": {
      "cylinders": [
        {"x": 75, "y": 25, "radius": 6.0},
        {"x": 75, "y": 45, "radius": 6.0},
        {"x": 105, "y": 25, "radius": 6.0},
        {"x": 105, "y": 45, "radius": 6.0}
      ]
    },
    "runtime": {
      "total_steps": 500,
      "warmup_steps": 100,
      "snapshot_every": 100
    },
    "outputs": {
      "output_profile": "analysis",
      "no_video": true,
      "paraview_output": false
    }
  },
  "factors": {
    "<factor1>": [<value1>, <value2>, ...],
    "<factor2>": [<value1>, <value2>, ...]
  }
}
```

**ルール:**
- `name` はタイムスタンプベースの一意な名前にする。**タイムスタンプは必ず `date +%Y%m%d_%H%M%S` を Bash で実行して動的に取得すること**
- `name` と `result_tag` はアンダースコアと英数字のみ使用（スラッシュ・スペース不可）
- `factors` のキーはテンプレートの有効なパラメータ名であること
- 総条件数（全ファクター値の直積）が 50 を超える場合は警告を表示し、条件数削減を提案すること
- ユーザーが `defaults` のパラメータ（ステップ数・ジオメトリ等）を指定した場合のみ上書きする

**指定可能なファクター例:**
```
particles:  particle_volume_fraction, particle_radius, density_ratio
physics:    rolling_friction, rolling_friction_coeff, sliding_friction,
            particle_attraction, attraction_strength, particle_repulsion, repulsion_strength
flow:       reynolds_number, pressure_drop
solver:     particle_fluid_coupling（"point_force" or "immersed_boundary"）
```

### ステップ 3: 確認表示

生成したスイープ設定を人間が読みやすい形で表示し、ユーザーに確認を求めてください。

```
## スイープ設定の確認

**スイープファイル:** configs/lbm_dem/sweeps/<name>.json
**総条件数:** N 条件（推定実行時間: 約 X 分）

### スイープファクター
| ファクター | 値 |
|-----------|-----|
| particle_volume_fraction | [0.1, 0.2, 0.3] |
| rolling_friction | [false, true] |

### 共通設定（defaults）
- ベースコンフィグ: fouling_supply.json
- ステップ数: 500
- ジオメトリ: 4本シリンダー

このスイープを実行しますか？
```

推定実行時間は「1条件あたり約2分 × 条件数」で概算してください。

### ステップ 4: ユーザー承認後の実行

ユーザーが承認したら、以下のコマンドを実行してください。スイープは長時間かかるため **タイムアウトを 7200 秒**に設定すること：

```bash
poetry run python src/runners/run_lbm_dem_sweep.py \
  --sweep configs/lbm_dem/sweeps/<name>.json
```

実行中の進捗出力（条件ごとの完了ログ）をそのまま表示してください。

### ステップ 5: 結果比較表の表示

実行完了後、スイープ結果ディレクトリ（`src/results/run_lbm_dem/<sweep_name>/`）の `status.tsv` を読み込み、以下の形式で比較表を表示してください：

```
## スイープ完了

**結果ディレクトリ:** src/results/run_lbm_dem/<sweep_name>/
**完了条件数:** N / N

### 条件別結果比較

| 条件 | ファウリング率 | 圧力降下(最終) | 通過粒子比 | 計算時間 |
|------|-------------|--------------|-----------|---------|
| pvf=0.1, rf=false | XX% | X.XXe-4 | X.XX | Xs |
| pvf=0.1, rf=true  | XX% | X.XXe-4 | X.XX | Xs |
| ...               | ... | ...      | ...  | ... |

### 傾向サマリー
- 最もファウリングが進んだ条件: ...
- 最も圧力降下が大きい条件: ...
```

`status.tsv` が存在しない場合は、各条件ディレクトリの `summary.json` を個別に読んで同等の表を作成してください。

## 注意事項

- 総条件数が多い（>20）場合は実行前に所要時間の見積もりを伝えてください
- 各呼び出しで常に新しいスイープファイルが生成されます（前回の修正機能はありません）
- `particle_fluid_coupling` の値は `"point_force"` または `"immersed_boundary"` のみ有効です
