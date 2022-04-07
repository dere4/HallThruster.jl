#try new potential solver, ie having extra cells on boundary
"""
    solve_potential!(; U::Matrix{Float64}, params::NamedTuple)
function to solve the potential equation derived from the generalized Ohm's law
and employing charge conservation using quasineutrality. Second derivatives approximated
with 2nd order central difference scheme, first derivatives with central difference.
Centers of computational mesh placed on edges of fluid mesh, therefore edges correspond to boundaries for fluid.
If required, interpolation is used to infer properties at mesh boundaries. Potential is a function of magnetic
field, anomalous and classical collision frequency, neutral and ion density as well as ion velocity, and electron density
and temperature leading to electron pressure. Solved by applying Thomas algorithm, which is of complexity O(n) and valid
if matrix tridiagonal and diagonally dominant. The latter assumption almost always holds unless there is a huge discontinuity
in either electron mobility or electron density.
"""
function solve_potential_edge!(U, params)
    #directly discretising the equation, conserves properties such as negative semidefinite etc...
    #add functionality for nonuniform cell size
    (;z_cell, index, ϕ_L, ϕ_R) = params
    nedges = length(z_cell) - 1

    (;pe, ne, μ, A, b, ϕ) = params.cache

    b[1] = ϕ_L
    b[end] = ϕ_R

    A.d[1] = 1.0
    A.du[1] = 0.0
    A.d[end] = 1.0
    A.dl[end] = 0.0

    mi = params.config.propellant.m

    @inbounds for i in 2:nedges-1

        ne⁻ = ne[i]
        ne⁺ = ne[i + 1]
        μ⁻ = μ[i]
        μ⁺ = μ[i+1]

        Δz = z_cell[i+1] - z_cell[i]

        Δz² = Δz^2

        ∇_neue = 0.0
        for Z in 1:params.config.ncharge
            ∇_neue += Z * (U[index.ρiui[Z], i + 1] - U[index.ρiui[Z], i]) / Δz / mi
        end

        #direct discretization, h/2 to each side
        A.dl[i - 1] = ne⁻ * μ⁻ / Δz²
        A.d[i] = -(ne⁻ * μ⁻ + ne⁺ * μ⁺) / Δz²
        A.du[i] = ne⁺ * μ⁺ / Δz²

        #source term, h/2 to each side
        b[i] = (μ⁻ * (pe[i - 1] + pe[i])/2 - (μ⁺ + μ⁻) * (pe[i] + pe[i + 1])/2 + μ⁺ * (pe[i + 1] + pe[i + 2])/2) / Δz² + ∇_neue

        #for order verification, change to simpler source term
        b[i] += params.config.source_potential(U, params, i)

    end
    return tridiagonal_solve!(ϕ, A, b)
end