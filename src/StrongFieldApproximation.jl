module StrongFieldApproximation

using ElectricFields
using Unitful
using UnitfulAtomic
using FastGaussQuadrature

using LinearAlgebra
using SparseArrays
using StaticArrays
using FillArrays

import DataStructures: Queue, enqueue!, dequeue!

using ProgressMeter
using TimerOutputs

include("threading.jl")
include("telescope_iterators.jl")

include("find_blocks.jl")

include("atom_models.jl")
include("momentum_grid.jl")

include("units.jl")
include("field.jl")
include("ion_propagators.jl")
include("volkov.jl")
include("windows.jl")

include("channels.jl")
include("system.jl")
include("diagrams.jl")
include("momenta.jl")

include("dyson_expansions.jl")
include("multi_channel_sfa.jl")

include("ionization_rates.jl")

end # module
