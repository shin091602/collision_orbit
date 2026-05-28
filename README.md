# collision_orbit

地球月系の円制限三体問題（CR3BP）における宇宙機の衝突軌道の数値計算。

月の重力特異点を正則化変換（LC変換・KS変換）によって除去し、
衝突軌道の構造をグリッドサーチによって体系的に調べる。

## 構成

| ファイル | 内容 |
|---|---|
| `src/lc_canonical_cr3bp.jl` | LC変換後の正準方程式（平面運動、状態変数 $u_1,u_2,w_1,w_2$） |
| `src/ks_canonical_cr3bp.jl` | KS変換後の正準方程式（空間運動、状態変数 $u_1\text{--}u_4,w_1\text{--}w_4$） |
| `src/utils_regularization.jl` | Cartesian座標 ↔ LC座標・KS座標の変換関数 |
| `scripts/jacobi_scan_collision_orbits/compute/run_planar.jl` | 平面LC系でヤコビ定数を走査し、衝突軌道を計算 |
| `scripts/jacobi_scan_collision_orbits/plot/plot_planar.jl` | 平面計算結果をCairoMakieで描画 |
| `scripts/jacobi_scan_collision_orbits/compute/run_spatial.jl` | 空間KS系でヤコビ定数を走査し、衝突軌道を計算 |
| `scripts/jacobi_scan_collision_orbits/plot/plot_spatial.jl` | 空間計算結果をCairoMakieで描画 |

生成データは `data/`、図は `results/` に保存する。どちらもGit管理対象外。

### 座標変換関数

```julia
# 平面（LC変換）
zeta = cart2lc(z, mu)   # z = [x, y, px, py] -> zeta = [u1, u2, w1, w2]
z    = lc2cart(zeta, mu)

# 空間（KS変換）
zeta = cart2ks(z, mu)   # z = [x, y, z, px, py, pz] -> zeta = [u1..u4, w1..w4]
z    = ks2cart(zeta, mu)
```

Cartesian座標はCR3BPの回転座標で、月の位置は $(1-\mu,0)$。
LC変換では月中心の相対位置

```math
x - (1-\mu) = u_1^2 - u_2^2,\qquad y = 2u_1u_2
```

を用いる。擬似時間は $s$ とし、平面LCでは $dt/ds = r_M = u_1^2+u_2^2$。

衝突点では $u=0$ となり、Cartesian運動量への逆変換は特異になる。
計算スクリプトでは衝突点のCartesian位置も保存し、特異な運動量成分は `NaN`
として残す。

## セットアップ

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 平面LC計算

```bash
# インタラクティブ
julia --project=. --threads=auto

# ヤコビ定数スキャン
julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_planar.jl

# バックグラウンド実行
tmux new-session -d -s calc "julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_planar.jl"
tmux attach -t calc
```

計算結果は次に保存される。

```text
data/jacobi_scan_planar.jld2
```

保存内容は `cart`, `lc`, `C_values`, `mu`。

## 空間KS計算

```bash
julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_spatial.jl
```

計算結果は次に保存される。

```text
data/jacobi_scan_spatial.jld2
```

保存内容は `cart`, `ks`, `C_values`, `mu`, `directions`, `N_ETA`, `N_PHI`, `N_PSI`。
衝突初期条件は KS 座標で `u = 0`, `|w|^2 = 8μ` とし、`w` の方向を
3-sphere 上の Hopf 座標グリッドでサンプリングする。

## プロット

```bash
julia --project=. scripts/jacobi_scan_collision_orbits/plot/plot_planar.jl
julia --project=. scripts/jacobi_scan_collision_orbits/plot/plot_spatial.jl
```

平面問題のため、Cartesian座標では `xy` 平面だけを出力する。
各ヤコビ定数ごとにフォルダを分け、同じフォルダ内では軌道ごとに色を変える。

```text
results/jacobi_scan/planar/C2.80/cart_xy.pdf
results/jacobi_scan/planar/C2.80/cart_xy_zoom.pdf
results/jacobi_scan/planar/C2.80/lc_u1.pdf
results/jacobi_scan/planar/C2.80/lc_u2.pdf
results/jacobi_scan/planar/C2.80/lc_w1.pdf
results/jacobi_scan/planar/C2.80/lc_w2.pdf
```

`cart_xy_zoom.pdf` は月近傍の

```text
0.8 < x < 1.2,  -0.2 < y < 0.2
```

を表示する。

空間問題では Cartesian 座標の `xyz` 3D表示、`xy`, `yz`, `zx` 射影、
それらを結合した通常版とzoom版、および各 KS 成分と結合版を出力する。

## 解析

```bash
julia --project=. scripts/jacobi_scan_collision_orbits/analyze/count_loops.jl
julia --project=. scripts/jacobi_scan_collision_orbits/analyze/plot_loop_counts.jl
julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/analyze/run_large_loop_scan.jl
```

`count_loops.jl` は月中心半径 `r_M <= 0.2` の領域内で、月中心 `xy` 射影の
巻き角が1周以上変化する軌道をループ軌道として数える。結果は
`data/jacobi_scan_loop_counts.jld2` と TSV に保存する。

`run_large_loop_scan.jl` は分析用の大規模計算を行い、`|u|^2 = 0.25` で
積分を停止する。既定では停止までの正則化座標と月中心位置も保存し、
結果を `data/jacobi_scan_loop_counts_large.jld2` に出力する。サンプル数は
環境変数 `N_ANGLES`, `N_ETA`, `N_PHI`, `N_PSI` で変更できる。軌道データを
保存しない場合は `SAVE_TRAJECTORIES=false` を指定する。

## テスト

現時点では専用の `test/` は未整備。新しく追加する場合はJulia標準の `Test`
を使い、以下を優先して確認する。

- `cart2lc` / `lc2cart` と `cart2ks` / `ks2cart` の往復精度
- 正則化ハミルトニアンの保存
- `run_planar.jl` の小規模スモークテスト
