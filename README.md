# collision_orbit

地球月系の円制限三体問題（CR3BP）における宇宙機の衝突軌道の数値計算。

KS変換（Kustaanheimo-Stiefel変換）で月の重力特異点を正則化し、
衝突軌道の構造をグリッドサーチによって体系的に調べる。

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
