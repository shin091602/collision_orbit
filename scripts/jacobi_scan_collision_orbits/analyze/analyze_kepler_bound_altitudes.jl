using JLD2
using Printf

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const EARTH_MOON_DISTANCE_KM = 384400.0
const MOON_RADIUS_KM = 1737.4
const ALTITUDES_KM = [100.0, 1000.0, 3000.0, 6000.0]
const RADII = (MOON_RADIUS_KM .+ ALTITUDES_KM) ./ EARTH_MOON_DISTANCE_KM

summary_names() = ["C"; ["bound_fraction_$(Int(a))km" for a in ALTITUDES_KM]]

function moon_radius(traj, mu; spatial)
    xM = traj[1, :] .- (1 - mu)
    y = traj[2, :]
    if spatial
        return sqrt.(xM .^ 2 .+ y .^ 2 .+ traj[3, :] .^ 2)
    end
    return sqrt.(xM .^ 2 .+ y .^ 2)
end

function interpolate_state(traj, r, target)
    i = findfirst(>=(target), r)
    isnothing(i) && return nothing
    i == 1 && return traj[:, 1]

    a = (target - r[i - 1]) / (r[i] - r[i - 1])
    return (1 - a) .* traj[:, i - 1] .+ a .* traj[:, i]
end

function kepler_energy(z, mu; spatial)
    if spatial
        v2 = z[4]^2 + (z[5] - (1 - mu))^2 + z[6]^2
    else
        v2 = z[3]^2 + (z[4] - (1 - mu))^2
    end
    return 0.5 * v2 - mu / moon_radius(reshape(z, :, 1), mu; spatial)[1]
end

function orbit_energies(traj, mu; spatial)
    r = moon_radius(traj, mu; spatial)
    energies = fill(NaN, length(RADII))

    for (j, target) in enumerate(RADII)
        z = interpolate_state(traj, r, target)
        isnothing(z) || (energies[j] = kepler_energy(z, mu; spatial))
    end

    return energies
end

function analyze(path; spatial)
    data = load(path)
    C_values = data["C_values"]
    mu = data["mu"]
    all_cart = data["cart"]

    energy = Dict{Float64, Matrix{Float64}}()
    bound = Dict{Float64, Matrix{Bool}}()
    summary = zeros(length(C_values), 1 + length(ALTITUDES_KM))

    for (i, C) in enumerate(C_values)
        trajs = all_cart[C]
        E = zeros(length(ALTITUDES_KM), length(trajs))

        for k in eachindex(trajs)
            E[:, k] = orbit_energies(trajs[k], mu; spatial)
        end

        B = E .< 0
        energy[C] = E
        bound[C] = B
        summary[i, :] .= [C; vec(sum(B, dims = 2) ./ length(trajs))]

        @printf("C = %.2f: bound fractions = %s\n", C, join(round.(summary[i, 2:end], digits = 3), ", "))
    end

    return energy, bound, summary, C_values, mu
end

function write_tsv(path, summary)
    open(path, "w") do io
        println(io, join(summary_names(), '\t'))
        for i in axes(summary, 1)
            println(io, join(summary[i, :], '\t'))
        end
    end
end

println("Planar")
planar_energy, planar_bound, planar_summary, planar_C, mu = analyze(
    joinpath(DATA_DIR, "large_jacobi_scan_planar.jld2");
    spatial = false,
)

println("Spatial")
spatial_energy, spatial_bound, spatial_summary, spatial_C, _ = analyze(
    joinpath(DATA_DIR, "large_jacobi_scan_spatial.jld2");
    spatial = true,
)

out = joinpath(DATA_DIR, "kepler_bound_altitudes.jld2")
save(out,
     "altitudes_km", ALTITUDES_KM,
     "radii", RADII,
     "summary_names", summary_names(),
     "planar_energy", planar_energy,
     "planar_bound", planar_bound,
     "planar_summary", planar_summary,
     "planar_C_values", planar_C,
     "spatial_energy", spatial_energy,
     "spatial_bound", spatial_bound,
     "spatial_summary", spatial_summary,
     "spatial_C_values", spatial_C,
     "mu", mu)

write_tsv(joinpath(DATA_DIR, "kepler_bound_altitudes_planar.tsv"), planar_summary)
write_tsv(joinpath(DATA_DIR, "kepler_bound_altitudes_spatial.tsv"), spatial_summary)

println("Saved to $(relpath(out, joinpath(@__DIR__, "../../..")))")
