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
