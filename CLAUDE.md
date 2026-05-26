# collision_orbit — CLAUDE.md

## プロジェクト概要

地球月系の円制限三体問題（CR3BP）における宇宙機の軌道計算。
KS変換（Kustaanheimo-Stiefel変換）で月の重力特異点を正則化し、
衝突軌道の構造をグリッドサーチによって体系的に調べる。

## 実行環境・ワークフロー

```
MacBook (VS Code)
    ↓ VS Code Remote SSH
Ubuntu (Intel Core i9)
    ├── src/       ← コア実装
    ├── scripts/   ← 実行スクリプト
    ├── results/   ← プロット・図（git管理外）
    ├── data/      ← 中間データ（git管理外）
    └── tmux       ← 長時間計算のセッション管理
```

- VS Code Remote SSH で MacBook から Ubuntu 上のコードを編集・実行
- 長時間のグリッドサーチ計算は tmux セッションで管理
- git は Ubuntu 側で管理し、MacBook からも参照

## 依存パッケージ（Project.toml）

| パッケージ | 用途 |
|---|---|
| DifferentialEquations | ODE/SDEソルバー（KS変換後の方程式を解く） |
| StaticArrays | 固定長配列（状態ベクトル高速化） |
| ForwardDiff | 自動微分（ヤコビアン等） |
| Roots | 方程式求解（断面条件など） |
| Interpolations | 補間 |
| JLD2 | 計算結果の保存・読み込み（HDF5互換） |
| Plots | プロット |
| Revise | 開発時のホットリロード |

## セットアップ

```bash
# Juliaパッケージのインストール
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 実行方法

### インタラクティブ実行

```bash
julia --project=. --threads=auto
```

### tmux でバックグラウンド実行（長時間計算）

```bash
# セッション作成
tmux new-session -d -s calc "julia --project=. --threads=auto main.jl"

# ログ確認
tmux attach -t calc

# セッション一覧
tmux ls
```

### スクリプト直接実行

```bash
julia --project=. --threads=auto script.jl
```

## コード規約

- 状態ベクトルは `StaticArrays` の `SVector` を使う（メモリ効率・速度）
- 計算結果は `JLD2` で `.jld2` ファイルに保存する
- グリッドサーチは並列化（`Threads.@threads` または `pmap`）を前提とする
- `Revise` は開発時のみ `using Revise` で読み込む

## ディレクトリ構成

```
collision_orbit/
├── src/          # コアライブラリ（CR3BP, KS変換, グリッドサーチ）
├── scripts/      # 実行スクリプト（計算ジョブ単位で1ファイル）
├── results/      # プロット図・可視化結果（.gitignore）
├── data/         # 中間データ・計算結果（.jld2, .h5, .csv）（.gitignore）
└── CLAUDE.md
```

## git 運用

- `src/` と `scripts/` のみgit管理
- `results/`・`data/` は `.gitignore` で除外
- MacBook の VS Code からリモートで編集し、Ubuntu 上で実行
