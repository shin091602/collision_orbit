# Repository Guidelines

## Project Structure & Module Organization

This repository studies collision orbits in the Earth-Moon CR3BP using LC and KS regularization. Core Julia routines live in `src/`: `lc_canonical_cr3bp.jl` and `ks_canonical_cr3bp.jl` define equations of motion, while `utils_regularization.jl` contains coordinate transforms. Run-oriented code is under `scripts/jacobi_scan_collision_orbits/`, split into `compute/` jobs and `plot/` scripts. Reference derivations belong in `ref/`. Generated datasets and figures go in ignored `data/` and `results/` directories.

## Build, Test, and Development Commands

- `julia --project=. -e 'using Pkg; Pkg.instantiate()'`: install dependencies from `Project.toml` and `Manifest.toml`.
- `julia --project=. --threads=auto`: start an interactive Julia session with the project environment and all available threads.
- `julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/compute/run_planar.jl`: run the planar Jacobi scan.
- `julia --project=. --threads=auto scripts/jacobi_scan_collision_orbits/plot/plot_planar.jl`: generate plots from saved scan data.
- `tmux new-session -d -s calc "julia --project=. --threads=auto <script.jl>"`: run long computations detached.

## Coding Style & Naming Conventions

Use idiomatic Julia with 4-space indentation. Prefer `snake_case` for functions and variables, `UPPER_CASE` for constants such as `W_COLL`, and descriptive physical names (`mu`, `C_values`, `tspan`). Use `StaticArrays` where fixed-size vectors improve performance. Use `joinpath(@__DIR__, ...)` for script-relative paths. Keep comments focused on equations, assumptions, and non-obvious numerical choices.

## Testing Guidelines

There is no committed `test/` harness yet. For new tests, use Julia's standard `Test` package under `test/runtests.jl`, and run them with `julia --project=. -e 'using Pkg; Pkg.test()'`. Prioritize transform round trips (`cart2lc`/`lc2cart`, `cart2ks`/`ks2cart`), conservation checks, and smoke tests for compute scripts. Name tests by behavior, for example `@testset "LC round trip"`.

## Commit & Pull Request Guidelines

Follow `RULES.md`: `<type>(<scope>): <summary>`, with an English imperative summary of 50 characters or less and no final period. Common types include `feat`, `fix`, `math`, `refactor`, `docs`, `test`, `sim`, `perf`, `chore`, and `wip`; common scopes include `lc`, `ks`, `cart`, `cr3bp`, `transform`, and `eom`. Use the body to explain why, and cite papers or equations in the footer when relevant. Pull requests should describe the numerical change, list commands run, and note generated files kept out of Git.

## Data & Configuration Tips

Do not commit large outputs, including `.jld2`, `.h5`, `.csv`, images, or movies. Keep reproducible parameters in scripts or documented notes, and write outputs beneath `data/` or `results/` using stable filenames.
