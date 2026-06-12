# 膜ファウリング研究ロードマップ

このドキュメントは、LBM-DEM連成モデルを使って解き明かすべき問いとその優先順位を整理したものです。

---

## モデルの強みと限界

**強み:**  
付着力(Hamaker引力/斥力)・転がり摩擦・体積分率・細孔形状を独立に振って、捕捉/通過の時系列を再現性よく比較できる。

**弱み(既知):**  
粒子は不透過境界として扱われておらず、ケーキ層の透過抵抗も未モデル化(→`docs/fouling_model/LBM_DEM_COUPLING_LIMITATIONS.md`)。

**照準:**  
閉塞の絶対値予測ではなく、**「ファウリングがどの様式で・いつ・なぜ始まるかの相対比較とメカニズム解明」**。

---

## Phase 1 — 今のスイープ完走後、追加計算ゼロで解く

`fouling_3d_20260609`(12ケース)の time_series.csv から読み取れる問い。

### 1. ファウリング様式の相図づくり ⭐ 本命

古典的 Hermia モデルでは、透過流束の減衰挙動を冪乗則でフィッティングし指数 n で様式を判別する:

| n | 様式 |
|---|---|
| 2 | 完全閉塞(Complete blocking) |
| 1 | 中間閉塞(Intermediate blocking) |
| 1.5 | 標準閉塞(Standard blocking) |
| 0 | ケーキろ過(Cake filtration) |

`removed_particles` の差分から透過流束プロキシ(`passed_rate = diff(removed) / diff(step)`)を計算し、各ケースの n を推定。**attraction_strength × particle_volume_fraction × cylinder_spec(gap)のパラメータ空間に様式マップを描く**。スイープ軸がそのまま相図の軸になる。

**使用データ:** `sweep_results/fouling_3d_20260609/<case>/analysis/time_series.csv`  
**実装:** `src/tools/analyze_lbm_dem_design_sweeps.py` に Hermia フィッティング関数を追加

### 2. 臨界付着条件(無次元数のスケーリング)

付着力と流体剥離力の競合を1つの無次元数で表す:

```
N_ad = Hamaker引力(attraction_strength) / (6πμ a U_max)
```

attraction 系列(a1–a3)の捕捉率を N_ad でプロットし、他のパラメータセットも同じ曲線に乗るかを検証。乗れば「臨界フラックス」の一般則が得られる。

**使用データ:** `summary.json`(captured_particles/total) + `metadata.json`(パラメータ)

### 3. 初期堆積の自己触媒効果

最初に円柱に付いた粒子が付着サイトを増やし後続の捕捉を加速するメカニズム。

- `particles_*.vtk` から「壁面直接付着」と「既付着粒子への付着」をステップごとに分類
- 捕捉確率が表面被覆率とともにどう上がるかを定量化(誘導期→加速期の遷移)

**使用データ:** `particles_*.vtk` または `particle_positions.npz`

---

## Phase 2 — 解析の追加で掘り下げる

スイープ結果を見てから優先度を再評価。

### 4. 堆積層の微細構造と付着力・摩擦の関係

rolling_friction(f1, f2)× attraction が堆積構造を「樹枝状(空隙大)」か「緻密」かに振り分ける。

- 粒子座標から**局所充填率・配位数・層表面粗さ**を計算
- 付着力/摩擦→層構造→流束低下速度の因果連鎖を定量化

**意義:** ケーキ抵抗未モデル化の弱点を「構造解析」で迂回する。

### 5. ブリッジング閉塞の統計

粒径/間隙比 ≳ 1/3 でアーチ閉塞が起きるとされる。ギャップ系列(g1, g2)と体積分率系列(c1, c2)の組み合わせで:

- 閉塞イベントを間隙ごとに検出
- 「何個の粒子が・どの確率でブリッジを形成するか」を統計化

---

## Phase 3 — モデル拡張で価値を大きく跳ねさせる

スイープ結果から得た物理的知見を土台にして実装判断。

### 6. 可逆 vs 不可逆ファウリング(逆洗シミュレーション)

ファウリング後に駆動力を反転 or 増速し、どの粒子が剥がれ・どれが残るかを観察。

- **実務直結:** 「物理洗浄で透水性はどこまで回復するか」に答えられる
- 付着力パラメータ×洗浄強度の**回復率マップ**
- ソルバー本体を変えず、流れ制御の切り替えだけで実装できる可能性が高い

### 7. ケーキ層透過抵抗(Brinkman 項の導入)

局所粒子充填率に応じた Darcy/Brinkman 抵抗を流体側に加算。

- `LBM_DEM_COUPLING_LIMITATIONS.md` で既に課題として記載
- 閉塞系の問い(1, 5)の定量精度が一段上がる
- Phase 2 の結果を見てから優先度を決める

---

## 進め方

```
今すぐ:  fouling_3d_20260609 (12ケース) を完走させる
         ↓
Phase 1: time_series + summary から Hermia 様式判別 + N_ad スケーリング解析
         ↓ 相図の空白を確認
追加計算: 不足パラメータ点を埋める次スイープを設計
         ↓
Phase 2: 粒子位置データから層構造・ブリッジング解析
         ↓ 知見が揃ったら
Phase 3: 逆洗 or Brinkman 拡張を実装
```

---

## 関連ファイル

| パス | 内容 |
|---|---|
| `docs/TOPOLOGY_OPTIMIZATION_PLAN.md` | サロゲート支援トポロジー最適化の研究計画(本ロードマップと指標を共有) |
| `docs/fouling_model/LBM_DEM_COUPLING_LIMITATIONS.md` | モデルの既知限界 |
| `docs/fouling_model/ALGORITHMS.md` | 数値手法の詳細 |
| `configs/lbm_dem/sweeps/manifests/fouling_3d_20260609.json` | 現行12ケースのマニフェスト |
| `src/tools/analyze_lbm_dem_design_sweeps.py` | スイープ解析スクリプト(Phase 1 拡張先) |
| `.claude/skills/analyze/` | `/analyze` スキル実装 |
