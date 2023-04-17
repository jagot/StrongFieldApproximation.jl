using StrongFieldApproximation
using ElectricFields
using Unitful
using UnitfulAtomic
using Test

# write your own tests here
@test 1 == 1
@testset "Telescope iterators" begin
    generate_reference_common(ref) =
        map(I -> [Tuple(I)...], findall(ref))

    generate_1d_reference(N, memory) = [[i] for i in 1:N]

    function generate_2d_reference(N, memory)
        ref = [0 < i - j ≤ memory
               for i in 1:N, j in 1:N]
        generate_reference_common(ref)
    end

    function generate_3d_reference(N, memory)
        ref = [i > j > k && 0 < i - j ≤ memory && 0 < j - k ≤ memory
               for i in 1:N, j in 1:N, k in 1:N]
        generate_reference_common(ref)
    end

    function generate_4d_reference(N, memory)
        ref = [i > j > k > l && 0 < i - j ≤ memory && 0 < j - k ≤ memory && 0 < k - l ≤ memory
               for i in 1:N, j in 1:N, k in 1:N, l in 1:N]
        generate_reference_common(ref)
    end

    function test_telescope_iterator(N, n, memory)
        v = 1:N
        ti = StrongFieldApproximation.TelescopeIterator(v, n, memory)
        # At the moment we don't know how to compute length(ti), so
        # collect(ti) does not work.
        data = Vector{Int}[]
        for e in ti
            push!(data, e)
        end
        ref = if n == 1
            generate_1d_reference(N, memory)
        elseif n == 2
            generate_2d_reference(N, memory)
        elseif n == 3
            generate_3d_reference(N, memory)
        elseif n == 4
            generate_4d_reference(N, memory)
        end

        @test data == ref
    end
    @testset "N = $N, n = $n, memory = $memory" for N = 1:7, n = 1:4, memory = [0,1,3,100]
        test_telescope_iterator(N, n, memory)
    end
end

@testset "Diagrams" begin
    function display_diagram(system, diagram, expected_momenta, expected_unique, expected_indeterminate)
        display(diagram)
        ions, unique_momenta, momenta, indeterminate_momenta, order = StrongFieldApproximation.analyze_diagram(system, diagram)

        expected_ions = first.(diagram.path)
        if first(diagram.path)[2] == 0 && length(diagram) > 1
            expected_ions = expected_ions[2:end]
        end

        @test ions == expected_ions
        @test momenta == expected_momenta
        @test unique_momenta == expected_unique
        @test indeterminate_momenta == expected_indeterminate
    end

    @field(F) do
        λ = 800.0u"nm"
        I₀ = 1e14u"W/cm^2"
        cycles = 4.0
        ϕ = π
        env = :cos²
    end

    ndt = 238

    Iₚ = 14u"eV" # "Krypton"

    # Elastic scattering off a Yukawa potential
    cc = StrongFieldApproximation.CoulombCoupling((𝐤,𝐩) -> yukawa_fourier(𝐩-𝐤, 1, 0, 1))
    dc = StrongFieldApproximation.DipoleCoupling(0.1, F)
    couplings=Matrix{StrongFieldApproximation.AbstractCoupling}[reshape([dc],1,1),reshape([cc],1,1)]

    ar = (F, ndt)
    channel = IonizationChannel(Iₚ, ar...)
    system = StrongFieldApproximation.System(repeat([channel], 10), nothing, couplings, ar...)

    for (path, expected_momenta, expected_unique, expected_indeterminate) in [
        ([(1,0)], [1], [(1,2)], []),
        ([(1,0),(1,0)], [1], [(1,2)], 1:1),
        ([(2,1),(1,0)], [1,1], [(1,3)], []),
        ([(2,2),(1,0)], [1,2], [(1,2),(2,3)], 2:2),
        ([(1,2),(2,2),(1,0)], [1,2,3], [(1,2),(2,3),(3,4)], 2:3),
        ([(1,0),(1,2),(2,2),(1,0)], [1,2,3], [(1,2),(2,3),(3,4)], 1:3),
        ([(1,1),(2,2),(1,0)], [1,1,2], [(1,3),(3,4)], 2:2),
        ([(1,0),(1,1),(2,2),(1,0)], [1,1,2], [(1,3),(3,4)], 1:2),
        ([(3,1),(1,1),(2,2),(1,0)], [1,1,1,2], [(1,4),(4,5)], 2:2),
        ([(3,0),(3,1),(1,1),(2,2),(1,0)], [1,1,1,2], [(1,4),(4,5)], 1:2),
        ([(1,2),(3,1),(1,1),(2,2),(1,0)], [1,2,2,2,3], [(1,2),(2,5),(5,6)], 2:3),
        ([(1,2),(1,2),(3,1),(2,2),(1,0)], [1,2,3,3,4], [(1,2),(2,3),(3,5),(5,6)], 2:4),
        ([(1,2),(1,2),(3,1),(1,1),(2,2),(1,0)], [1,2,3,3,3,4], [(1,2),(2,3),(3,6),(6,7)], 2:4),
        ([(1,2),(10,2),(3,2),(1,2),(2,2),(1,0)], [1,2,3,4,5,6], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7)], 2:6),
        ([(1,0),(1,2),(10,2),(3,2),(1,2),(2,2),(1,0)], [1,2,3,4,5,6], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7)], 1:6),
        ([(1,0),(1,1),(10,1),(3,1),(1,1),(2,1),(1,0)], [1,1,1,1,1,1], [(1,7)], 1:1)
    ]
        diagram = Diagram(path, system)
        display_diagram(system, diagram, expected_momenta, expected_unique, expected_indeterminate)
    end
end
