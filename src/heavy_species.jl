
function update_heavy_species!(dU, U, params, t)
    ####################################################################
    #extract some useful stuff from params

    (;index, z_edge) = params
    (;
        source_neutrals, source_ion_continuity, source_ion_momentum,
        electron_pressure_coupled, propellant, scheme
    ) = params.config

    (;ue, μ, F, UL, UR) = params.cache

    ncells = size(U, 2)

    mi = propellant.m

    ##############################################################
    #FLUID MODULE

    #fluid BCs now in update_values struct

    ncharge = params.config.ncharge

    compute_fluxes!(F, UL, UR, U, params)

    @inbounds for i in 2:ncells-1
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        dU[index.ρn, i] = (F[index.ρn, left] - F[index.ρn, right]) / Δz

        # User-provided neutral source term
        dU[index.ρn, i] += source_neutrals(U, params, i)

        for Z in 1:ncharge

            dU[index.ρi[Z]  , i] = (F[index.ρi[Z],   left] - F[index.ρi[Z],   right]) / Δz
            dU[index.ρiui[Z], i] = (F[index.ρiui[Z], left] - F[index.ρiui[Z], right]) / Δz
            # User-provided source terms
            dU[index.ρi[Z],   i] += source_ion_continuity[Z](U, params, i)
            dU[index.ρiui[Z], i] += source_ion_momentum[Z  ](U, params, i)
        end

        apply_ion_acceleration!(dU, U, params, i)
        apply_reactions!(dU, U, params, i)

    end

    return nothing
end