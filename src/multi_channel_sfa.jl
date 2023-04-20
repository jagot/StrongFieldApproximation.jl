function ionization!(v, system, 𝐩, 𝐀, i, tmp)
    p = kinematic_momentum(𝐩, 𝐀[i])
    for j in eachindex(v)
        tmp[j] = source_term(system.ionization_channels[j], i, p)
    end
    mul!(v, ion_mapping(system.ions, i), tmp)
end

macro swap!(x,y,tmp)
    quote
        $(esc(tmp)) = $(esc(x))
        $(esc(x)) = $(esc(y))
        $(esc(y)) = $(esc(tmp))
    end
end

function multi_channel_sfa!(amplitudes, system, diagram, iref, 𝐤=nothing;
                            window=flat_window(), imin=1, verbosity=0,
                            to=TimerOutput())
    _, unique_momenta, momenta, indeterminate_momenta, order = @timeit to "Analyze diagram" analyze_diagram(system, diagram)
    @timeit to "Allocations" begin
        weight = (-im*system.dt)^order
        if verbosity > 0
            println()
            @info "Integrating diagram up to" iref system diagram unique_momenta momenta indeterminate_momenta order weight 𝐤
        end

        𝐩s = complex(zeros(momentum_type(system, 𝐤), length(unique_momenta)))
        prefactors = ones(complex(eltype(system.t)), length(unique_momenta))
        if !isnothing(𝐤)
            𝐩s[1] = 𝐤
        end
        verbosity > 2 && @show 𝐩s

        𝐀 = system.𝐀
        amplitudes .= false
        u = similar(amplitudes)
        v = similar(amplitudes)
        tmp = similar(amplitudes)
        ctT = complex(eltype(system.t))

        is = vcat(iref, zeros(Int, order))

        memory = length(window)
        if memory < iref
            window = vcat(window, zeros(eltype(window), iref))
        end
    end

    @timeit to "Time loop" begin
        for i in TelescopeIterator(max(1,imin):iref-1, order, memory)
            windw = @timeit to "Window" prod(window[(i[1] + 1) .- i])
            iszero(windw) && continue

            is[2:end] .= i

            u .= false
            v .= false

            @timeit to "Evaluate momenta" evaluate_momenta!(𝐩s, prefactors, system, unique_momenta, indeterminate_momenta, is)
            preprod = prod(prefactors)

            @timeit to "Ionization" ionization!(u, system, 𝐩s[end], 𝐀, is[end], v)

            @timeit to "Interior orders" begin
                # We have to perform this loop in chronological order, because
                # of the pesky time-ordering operator. Since the Goldstone
                # diagrams are stored with their vertices in
                # anti-chronological order, we thus have to loop backwards.
                for j = order:-1:1
                    a,b = is[j], is[j+1]

                    # First, if there is any interaction (beyond ionization)
                    # at this vertex, we affect this by simply acting with it
                    # on the current amplitudes.

                    which = diagram.path[j][2]
                    if which ≠ 0
                        @timeit to "Interaction" begin
                            𝐀ᵢ = 𝐀[b]
                            𝐤ᵢ = kinematic_momentum(𝐩s[momenta[j]], 𝐀ᵢ)
                            𝐩ᵢ = kinematic_momentum(𝐩s[momenta[j+1]], 𝐀ᵢ)

                            apply_interaction!(v, system.ions, system.couplings[which],
                                u, 𝐤ᵢ, 𝐩ᵢ, b)
                            @timeit to "Swap solutions" begin
                                @swap!(u,v,tmp)
                            end
                        end
                    end

                    # In-between interactions, it is assumed that the
                    # interaction-free ion+photoelectron propagator
                    # factorizes:
                    #
                    # U⁽⁰⁾(a,b) = Uᵢₒₙ(a,b) Uₑₗ(a,b)
                    #
                    # and we may thus propagate the electrons independently,
                    # and then mix their channel amplitudes by the "rotation"
                    # accumulated over the time interval by propagating the
                    # ions (including external field).

                    # Propagate electrons from b → a
                    @timeit to "Volkov propagation" begin
                        Sₑₗ = @timeit to "Volkov phase" volkov_phase(𝐩s[j], system.volkov, a, b)
                        eSₑₗ = exp(-im*Sₑₗ)
                        u .*= eSₑₗ
                    end

                    @timeit to "Ion/channel propagation" begin
                        # Propagate ions from b → a
                        u .*= ion_propagation(system.ions, :, a, b)
                    end
                end
            end

            order > 1 && last(first(diagram)) == 0 && error("Recombination not yet implemented in multi-channel case")

            amplitudes .+= u .* (windw*preprod)
        end
    end

    @timeit to "Weighting" lmul!(weight, amplitudes)
end

function multi_channel_photoelectron_spectrum(k::AbstractArray{T},
                                              system::System, diagram::Diagram;
                                              iref=length(system.t),
                                              verbosity=1, kwargs...) where T
    verbosity == 1 && @info "Ostensibly multi-channel photoelectrum spectrum calculation" system diagram length(k)

    nc = num_channels(system)
    sk = size(k)

    cT = complex(eltype(T))
    c = zeros(cT, (nc, sk...))
    to = TimerOutput()
    p = Progress(length(k))
    @timeit to "k loop" begin
        threaded_range_loop(CartesianIndices(k)) do I
            tok = TimerOutput()
            multi_channel_sfa!(view(c, :, I), system, diagram, iref, k[I]; verbosity=verbosity-10, to=tok, kwargs...)
            ProgressMeter.next!(p)
            merge!(to, tok, tree_point=["k loop"])
        end
    end
    TimerOutputs.complement!(to)
    0 < verbosity < 4 && print_timer(to)
    c
end

export multi_channel_photoelectron_spectrum
