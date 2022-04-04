
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

    if scheme.WENO #with WENO 5 only going up to 2nd last to boundary cells
        @inbounds for i in 3:ncells-3 #only fluxes ncells - 1
            # Handle neutrals
            F[index.ρn, i] = WENO5_compute_fluxes(
                F[index.ρn, i-2], F[index.ρn, i-1], F[index.ρn, i], F[index.ρn, i+1], F[index.ρn, i+2])
            #Handle ions
            for Z in 1:ncharge
                @views @. F[index.ρi[Z]:index.ρiui[Z], i] = WENO5_compute_fluxes(
                    SA[F[index.ρi[Z], i-2], F[index.ρiui[Z], i-2]],
                    SA[F[index.ρi[Z], i-1], F[index.ρiui[Z], i-1]],
                    SA[F[index.ρi[Z], i],   F[index.ρiui[Z], i]],
                    SA[F[index.ρi[Z], i+1], F[index.ρiui[Z], i+1]],
                    SA[F[index.ρi[Z], i+2], F[index.ρiui[Z], i+2]],
                )
                
            end
        end

        
        #for ion OVS test do periodic boundaries effectively
        #cell 1 #################
        F[index.ρn, 1] = WENO5_compute_fluxes(
                F[index.ρn, end-2], F[index.ρn, end-1], F[index.ρn, 1], F[index.ρn, 2], F[index.ρn, 3])
        #Handle ions
        for Z in 1:ncharge
            @views @. F[index.ρi[Z]:index.ρiui[Z], 1] = WENO5_compute_fluxes(
                SA[F[index.ρi[Z], end-2], F[index.ρiui[Z], end-2]],
                SA[F[index.ρi[Z], end-1], F[index.ρiui[Z], end-1]],
                SA[F[index.ρi[Z], 1],   F[index.ρiui[Z], 1]],
                SA[F[index.ρi[Z], 2], F[index.ρiui[Z], 2]],
                SA[F[index.ρi[Z], 3], F[index.ρiui[Z], 3]],
            )
            
        end
        #cell 2 #################
        F[index.ρn, 1] = WENO5_compute_fluxes(
                F[index.ρn, end-1], F[index.ρn, 1], F[index.ρn, 2], F[index.ρn, 3], F[index.ρn, 4])
        #Handle ions
        for Z in 1:ncharge
            @views @. F[index.ρi[Z]:index.ρiui[Z], 1] = WENO5_compute_fluxes(
                SA[F[index.ρi[Z], end-1], F[index.ρiui[Z], end-1]],
                SA[F[index.ρi[Z], 1], F[index.ρiui[Z], 1]],
                SA[F[index.ρi[Z], 2], F[index.ρiui[Z], 2]],
                SA[F[index.ρi[Z], 3], F[index.ρiui[Z], 3]],
                SA[F[index.ρi[Z], 4], F[index.ρiui[Z], 4]],
            )
            
        end
        #cell end-2 #################
        F[index.ρn, 1] = WENO5_compute_fluxes(
                F[index.ρn, end-1], F[index.ρn, 1], F[index.ρn, 2], F[index.ρn, 3], F[index.ρn, 4])
        #Handle ions
        for Z in 1:ncharge
            @views @. F[index.ρi[Z]:index.ρiui[Z], 1] = WENO5_compute_fluxes(
                SA[F[index.ρi[Z], end-1], F[index.ρiui[Z], end-1]],
                SA[F[index.ρi[Z], 1], F[index.ρiui[Z], 1]],
                SA[F[index.ρi[Z], 2], F[index.ρiui[Z], 2]],
                SA[F[index.ρi[Z], 3], F[index.ρiui[Z], 3]],
                SA[F[index.ρi[Z], 4], F[index.ρiui[Z], 4]],
            )
            
        end


    end

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