# collision_orbit

地球月系の円制限三体問題（CR3BP）における宇宙機の衝突軌道の数値計算。

月の重力特異点を正則化変換（LC変換・KS変換）によって除去し、
衝突軌道の構造をグリッドサーチによって体系的に調べる。

## ライブラリ構成（`src/`）

| ファイル | 内容 |
|---|---|
| `lc_canonical_cr3bp.jl` | LC変換後の正準方程式（平面運動、状態変数 $\bar{x}_1,\bar{x}_2,\bar{p}_1,\bar{p}_2$） |
| `ks_canonical_cr3bp.jl` | KS変換後の正準方程式（空間運動、状態変数 $u_1\text{–}u_4,\,w_1\text{–}w_4$） |
| `utils_regularization.jl` | 直交座標 ↔ LC座標・KS座標の変換関数 |

### 座標変換関数

```julia
# 平面（LC変換）
zeta = cart2lc(z, mu)   # z = [q1, q2, p1, p2]  →  zeta = [ū1, ū2, p̄1, p̄2]
z    = lc2cart(zeta, mu)

# 空間（KS変換）
zeta = cart2ks(z, mu)   # z = [q1, q2, q3, p1, p2, p3]  →  zeta = [u1..u4, w1..w4]
z    = ks2cart(zeta, mu)
```

座標系はすべて月中心の回転座標系（擬似時間 $s$、$dt/ds = r$）。
正準方程式の導出は [`ref/canonical_equation_LC_KS.md`](ref/canonical_equation_LC_KS.md) を参照。

## セットアップ

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 実行

```bash
# インタラクティブ
julia --project=. --threads=auto

# バックグラウンド（tmux）
tmux new-session -d -s calc "julia --project=. --threads=auto main.jl"
tmux attach -t calc
```
