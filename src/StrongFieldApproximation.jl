module StrongFieldApproximation

using ElectricFields
using Unitful
using UnitfulAtomic

using FastGaussQuadrature
using IntervalSets

using LinearAlgebra
using StaticArrays

using ProgressLogging

include("threading.jl")
include("telescope_iterators.jl")

include("atom_models.jl")
include("momentum_grid.jl")

include("units.jl")
include("field.jl")
include("automatic_integration.jl")
include("volkov.jl")

include("channels.jl")
include("dyson_expansions.jl")

include("ionization_rates.jl")

end # module
