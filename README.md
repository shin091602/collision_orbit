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
| `scripts/jacobi_scan_collision_orbits/compute/run_planar.jl` | 平面LC系でヤコビ定数を走査し、衝突軌道を計算（小規模） |
| `scripts/jacobi_scan_collision_orbits/compute/run_large_planar.jl` | 平面LC系でヤコビ定数を走査し、衝突軌道を計算（大規模・解析用） |
| `scripts/jacobi_scan_collision_orbits/compute/run_spatial.jl` | 空間KS系でヤコビ定数を走査し、衝突軌道を計算（小規模） |
| `scripts/jacobi_scan_collision_orbits/compute/run_large_spatial.jl` | 空間KS系でヤコビ定数を走査し、衝突軌道を計算（大規模・解析用） |
| `scripts/jacobi_scan_collision_orbits/plot/plot_planar.jl` | 平面計算結果をCairoMakieで描画 |
| `scripts/jacobi_scan_collision_orbits/plot/plot_spatial.jl` | 空間計算結果をCairoMakieで描画 |
| `scripts/jacobi_scan_collision_orbits/analyze/analyze_loop_classification.jl` | ループ分類指標を計算・保存 |
| `scripts/jacobi_scan_collision_orbits/analyze/plot_loop_classification.jl` | ループ分類結果をプロット |
| `scripts/jacobi_scan_collision_orbits/analyze/analyze_kepler_bound_altitudes.jl` | 高度別ケプラーエネルギー拘束割合を計算・保存 |
| `scripts/jacobi_scan_collision_orbits/analyze/plot_kepler_bound_altitudes.jl` | ケプラーエネルギー拘束割合をプロット |

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

# ヤコビ定数スキャン（小規模）
julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_planar.jl

# 解析用大規模スキャン（C = 2.90:0.01:3.20、2000方向）
julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_large_planar.jl

# バックグラウンド実行
tmux new-session -d -s calc "julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_large_planar.jl"
tmux attach -t calc
```

### 衝突初期条件（平面）

正則化ハミルトニアンの衝突条件は $u=0$、$|w|^2 = 8\mu$。
初期速度方向を $\alpha \in [0, 2\pi)$ で一様にサンプリングする：

```math
w_1 = \sqrt{8\mu}\cos\alpha,\quad w_2 = \sqrt{8\mu}\sin\alpha
```

積分は $|u|^2 \geq 0.22$（月から十分離れた距離）で打ち切る。

計算結果は `data/jacobi_scan/large_jacobi_scan_planar.jld2` に保存される。
保存内容は `cart`, `lc`, `C_values`, `mu`。

## 空間KS計算

```bash
# 小規模
julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_spatial.jl

# 解析用大規模スキャン（C = 2.90:0.01:3.20、N_ETA=10, N_PHI=44, N_PSI=44）
julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_large_spatial.jl
```

### 衝突初期条件（空間）

KS座標の衝突条件は $u=0$、$|w|^2 = 8\mu$（LC条件の4次元類比）。
速度方向 $w/|w|$ を3次元球面 $S^3$ 上の Hopf 座標グリッドで一様サンプリングする：

```math
\rho_1 = \sqrt{(i-0.5)/N_\eta},\quad \rho_2 = \sqrt{1-\rho_1^2}
```
```math
w = |w|\,[\rho_1\cos\phi,\ \rho_1\sin\phi,\ \rho_2\cos\psi,\ \rho_2\sin\psi]
```

$\rho_1^2$ を $[0,1]$ 上で一様にとることで $S^3$ 上の面積素を均等にカバーする。

計算結果は `data/jacobi_scan/large_jacobi_scan_spatial.jld2` に保存される。
保存内容は `cart`, `ks`, `C_values`, `mu`, `directions`, `N_ETA`, `N_PHI`, `N_PSI`。

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
...
```

空間問題では Cartesian 座標の `xyz` 3D表示、`xy`, `yz`, `zx` 射影、
それらを結合した通常版とzoom版、および各 KS 成分と結合版を出力する。

## 解析

解析スクリプトはいずれも大規模計算結果（`data/jacobi_scan/large_jacobi_scan_*.jld2`）を入力とする。
`analyze_*.jl` が JLD2 と TSV の集計データを作成し、対応する `plot_*.jl` が
`results/jacobi_scan/analysis/` 以下に要約図を出力する。

### ループ分類

```bash
julia --project=. scripts/jacobi_scan_collision_orbits/analyze/analyze_loop_classification.jl
julia --project=. scripts/jacobi_scan_collision_orbits/analyze/plot_loop_classification.jl
```

各衝突軌道が月近傍でどのような形状をとるかを分類する。
分析対象区間は月中心距離 $r_M \in [0.02,\ 0.2]$（無次元単位、地球月距離=1）。
軌道が $r_M = 0.02$ に到達しない場合は未計測として扱い、`n_measured` に含めない。

#### 方位角ループ（azimuthal loop）

月中心からの方位角 $\theta = \mathrm{atan2}(y_M, x_M)$ を unwrap し、

```math
\Delta\theta_{\max} = \max_t |\theta(t) - \theta(t_0)|
```

とするとき、$\Delta\theta_{\max} / 2\pi \geq 1$ ならば方位角ループと判定する。
月の周りを少なくとも1周以上まわる軌道に対応する。

#### 半径ループ（radial loop）

分析区間内での地球中心距離 $r_E$ の局所的極大の個数が2以上のとき、
半径ループと判定する。
軌道が月近傍で地球方向に対して振動（行き来）を繰り返すことを示す。

#### 結合ループ（combined loop）

方位角ループまたは半径ループのどちらか一方でも満たせば combined loop とする。

#### 出力

```text
data/jacobi_scan/loop_classification.jld2
data/jacobi_scan/loop_classification_planar.tsv
data/jacobi_scan/loop_classification_spatial.tsv
results/jacobi_scan/analysis/loop_classification_summary.pdf
```

TSV の列は `C`, `n_orbits`, `n_measured`, `azimuthal_fraction`, `radial_fraction`, `combined_fraction`。
各 `*_fraction` は全軌道数 `n_orbits` を分母にした割合である。

### ケプラーエネルギー拘束高度解析

```bash
julia --project=. scripts/jacobi_scan_collision_orbits/analyze/analyze_kepler_bound_altitudes.jl
julia --project=. scripts/jacobi_scan_collision_orbits/analyze/plot_kepler_bound_altitudes.jl
```

衝突軌道が特定の高度を通過した時点での「月に対するケプラーエネルギー」を計算し、
負エネルギー（月重力に捕捉された状態）の軌道割合をヤコビ定数の関数として調べる。

#### ケプラーエネルギーの定義

月中心距離 $r_M$ が対象高度に達した時刻での状態ベクトルを補間し、

```math
E_K = \frac{1}{2}|\mathbf{v}_{\rm rel}|^2 - \frac{\mu}{r_M}
```

を計算する。ここで $\mathbf{v}_{\rm rel}$ は**慣性系**における月に対する相対速度である。
回転座標での速度成分を平面では $(v_x, v_y)$、空間では $(v_x, v_y, v_z)$ とする。
月の慣性系速度は $(0,\ 1-\mu, 0)$（回転角速度 $\omega=1$）なので、

```math
|\mathbf{v}_{\rm rel}|^2 = v_x^2 + (v_y - (1-\mu))^2
```

空間KSの場合はこれに $v_z^2$ を加える。

$E_K < 0$ であれば、その高度で軌道が月の二体問題的重力に捕捉されている
（楕円軌道的状態）ことを意味する。

#### 解析高度

月面からの高度 100, 1000, 3000, 6000 km の4点。
地球月距離 $384400\ \rm km$ で正規化した無次元半径に換算して使用する。

```math
r = \frac{R_{\rm Moon} + h_{\rm alt}}{d_{EM}} = \frac{1737.4 + h}{384400}
```

#### 出力

```text
data/jacobi_scan/kepler_bound_altitudes.jld2
data/jacobi_scan/kepler_bound_altitudes_planar.tsv
data/jacobi_scan/kepler_bound_altitudes_spatial.tsv
results/jacobi_scan/analysis/kepler_bound_altitudes_summary.pdf
```

TSV の列は `C`, `bound_fraction_100km`, `bound_fraction_1000km`,
`bound_fraction_3000km`, `bound_fraction_6000km`。
各割合は対象高度に到達し、かつ $E_K < 0$ となった軌道数を全軌道数で割った値である。

## テスト

現時点では専用の `test/` は未整備。新しく追加する場合はJulia標準の `Test`
を使い、以下を優先して確認する。

- `cart2lc` / `lc2cart` と `cart2ks` / `ks2cart` の往復精度
- 正則化ハミルトニアンの保存
- `run_planar.jl` の小規模スモークテスト
