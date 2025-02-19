module HallThruster

using StaticArrays
using CSV
using DataFrames
using OrdinaryDiffEq
using DiffEqBase
using LoopVectorization
using LinearAlgebra
using FileIO
using DiffEqCallbacks
using SparseArrays
using PartialFunctions
using QuadGK
using DelimitedFiles: readdlm, writedlm
using DocStringExtensions
using Unitful
using SpecialFunctions

# Packages used for making plots
using Measures: mm
using RecipesBase

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")
const LANDMARK_FOLDER = joinpath(PACKAGE_ROOT, "landmark")
const LANDMARK_RATES_FILE = joinpath(LANDMARK_FOLDER, "landmark_rates.csv")

include("utilities/utility_functions.jl")
include("utilities/transition_functions.jl")

include("physics/physicalconstants.jl")
include("physics/gas.jl")
include("physics/conservationlaws.jl")
include("physics/fluid.jl")
include("physics/thermodynamics.jl")
include("physics/walls.jl")
include("physics/electrontransport.jl")

include("numerics/finite_differences.jl")
include("numerics/limiters.jl")
include("numerics/flux.jl")

include("collisions/reactions.jl")
include("collisions/ionization.jl")
include("collisions/excitation.jl")
include("collisions/elastic.jl")
include("collisions/collision_frequencies.jl")

include("simulation/initialization.jl")
include("simulation/geometry.jl")
include("simulation/postprocess.jl")
include("simulation/boundaryconditions.jl")
include("simulation/potential.jl")
include("simulation/heavy_species.jl")
include("simulation/electronenergy.jl")
include("simulation/sourceterms.jl")
include("simulation/update_values.jl")
include("simulation/configuration.jl")
include("simulation/restart.jl")
include("simulation/simulation.jl")
include("visualization/plotting.jl")
include("visualization/recipes.jl")

export time_average, Xenon, Krypton

end # module
