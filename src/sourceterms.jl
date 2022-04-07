function apply_reactions!(dU::AbstractArray{T}, U::AbstractArray{T}, params, i::Int64) where T
    (;index, reactions, index, species_range_dict) = params

    ne = electron_density(U, params, i)
    ϵ  = U[index.nϵ, i] / ne

    @inbounds for rxn in reactions
        reactant_index = species_range_dict[rxn.reactant.symbol][1]
        product_index  = species_range_dict[rxn.product.symbol ][1]
        ρ_reactant = U[reactant_index, i]
        ρdot = reaction_rate(rxn, ne, ρ_reactant, ϵ)
        dU[reactant_index, i] -= ρdot
        dU[product_index, i]  += ρdot
    end
end

@inline reaction_rate(rxn, ne, n_reactant, ϵ) = rxn.rate_coeff(ϵ) * n_reactant * ne

function apply_ion_acceleration!(dU, U, params, i)
    index = params.index
    (;∇ϕ, ue, μ) = params.cache
    coupled = params.config.electron_pressure_coupled
    mi = params.config.propellant.m

    @inbounds for Z in 1:params.config.ncharge
        Q_accel = coupled * ue[i] / μ[i] + (1 - coupled) * ∇ϕ[i]
        Q_accel = -Z * e * U[index.ρi[Z], i] / mi * Q_accel
        dU[index.ρiui[Z], i] += Q_accel
    end
end

function source_electron_energy!(Q, U, params, i)
    Q[params.index.nϵ] = source_electron_energy(U, params, i)
end

function source_electron_energy(U, params, i)
    index = params.index

    mi = params.config.propellant.m
    ne = params.cache.ne[i]
    ue = params.cache.ue[i]
    ∇ϕ = params.cache.∇ϕ[i]

    nn = U[index.ρn, i] / mi
    K = params.config.collisional_loss_model(U, params, i)
    W = params.config.wall_loss_model(U, params, i)

    ohmic_heating      = ne * ue * ∇ϕ
    wall_losses        = ne * W
    collisional_losses = ne * nn * K

    return ohmic_heating - wall_losses - collisional_losses
end