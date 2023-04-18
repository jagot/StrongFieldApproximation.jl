using StrongFieldApproximation
using ElectricFields
using Unitful
using UnitfulAtomic
using DelimitedFiles
using Test

reffile(name) = joinpath(dirname(@__FILE__), "references", name)

function test_approx_eq(a, b; on_fail::Union{Nothing,Function}=nothing, kwargs...)
    size(a) == size(b) || throw(DimensionMismatch("Cannot compare objects of sizes $(size(a)) and $(size(b))"))
    if !isapprox(a, b; kwargs...)
        @error "Approximate equality failed:"
        na = norm(a)
        nb = norm(b)
        Δ = norm(a-b)
        relΔ = Δ/max(na,nb)

        println("   |a|", na)
        println("   |b|", nb)
        println("Abs. Δ", Δ)
        println("Rel. Δ", relΔ)

        isnothing(on_fail) || on_fail()
    end

    @test isapprox(a, b; kwargs...)
end

@testset "StrongFieldApproximation.jl" begin
    include("telescope_iterators.jl")
    include("diagrams.jl")
    include("rescattered_ati.jl")
end
