using LinearAlgebra

function lc_canonical_cr3bp(dzeta, zeta, params, t)
    mu = params[1]
    C = params[2]
    u = zeta[1:2]
    w = zeta[3:4]

    r = u[1]^2 + u[2]^2
    r_E = sqrt((u[1]^2 - u[2]^2 + 1)^2 + 4 * u[1]^2 * u[2]^2)

    dzeta[1] = w[1]/4 + (r - 1 + mu) * u[2] / 2
    dzeta[2] = w[2]/4 - (r + 1 - mu) * u[1] / 2

    dzeta[3] = (-C * u[1]
                - w[1] * u[1] * u[2]
                + w[2] * (3 * u[1]^2 + u[2]^2 + 1 - mu) / 2
                + 2 * u[1] * (1 - mu) / r_E
                - 2 * u[1] * r * (r + 1) * (1 - mu) / r_E^3)
    dzeta[4] = (-C * u[2]
                - w[1] * (u[1]^2 + 3 * u[2]^2 - 1 + mu) / 2
                + w[2] * u[1] * u[2]
                + 2 * u[2] * (1 - mu) / r_E
                + 2 * u[2] * r * (1 - r) * (1 - mu) / r_E^3) 
end
