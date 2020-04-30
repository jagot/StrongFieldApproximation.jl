module StrongFieldApproximation

using ElectricFields
using Unitful
using UnitfulAtomic
using FastGaussQuadrature

include("atom_models.jl")
include("units.jl")
include("field.jl")
include("trapz.jl")
include("propagate.jl")

end # module
