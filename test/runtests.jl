using SFA
using Base.Test

# write your own tests here
@test 1 == 1
include("test_trapz.jl")
include("test_sfa.jl")
