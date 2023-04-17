# * Integrate diagrams

function ionization(system::System, diagram::Diagram, 𝐩, 𝐀, i)
    α,which = diagram.path[end]
    @assert which == 0
    source_term(system.ionization_channels[α],
                i,
                kinematic_momentum(𝐩, 𝐀[i]))
end

function recombination(system::System, diagram::Diagram, 𝐩, 𝐀, i)
    α,which = first(diagram)
    if which == 0 && length(diagram) > 1
        d = system.ionization_channels[α].st.d
        conj(d(kinematic_momentum(𝐩, 𝐀[i])))
    else
        true
    end
end

function integrate_diagram(::Type{Amp}, system::System, diagram::Diagram, iref, 𝐤=nothing; memory=typemax(Int), imin=1,
                           to=TimerOutput(), verbosity=1) where Amp
    ions, unique_momenta, momenta, indeterminate_momenta, order = @timeit to "Analyze diagram" analyze_diagram(system, diagram)

    @timeit to "Allocations" begin
        weight = (-im*system.dt)^order

        verbosity > 1 && @info "Integrating diagram up to" iref system diagram ions unique_momenta momenta indeterminate_momenta order weight 𝐤

        Eᵢₒₙₛ = [system.ionization_channels[ion].E for ion in ions]
        𝐩s = complex(zeros(momentum_type(system, 𝐤), length(unique_momenta)))
        prefactors = ones(complex(eltype(system.t)), length(unique_momenta))
        if !isnothing(𝐤)
            𝐩s[1] = 𝐤
        end

        𝐀 = system.𝐀
        amplitude = complex(zero(Amp))
        ctT = complex(eltype(system.t))

        is = vcat(iref, zeros(Int, order))
    end

    @timeit to "Time loop" begin
        for i in TelescopeIterator(max(1,imin):iref-1, order, memory)
            is[2:end] .= i
            # is = vcat(iref, i)
            # is = (iref,i...)
            # for i in 1:iref-1, j in 1:i-1
            #     is = (iref,i,j)
            # println(is)

            @timeit to "Evaluate momenta" evaluate_momenta!(𝐩s, prefactors, system, unique_momenta, indeterminate_momenta, is)
            verbosity > 10 && @show 𝐩s prefactors

            @timeit to "Evaluate propagators" begin
                Sᵢₒₙ = zero(ctT)
                Sₑₗ = zero(ctT)
                for j = 1:order
                    a,b = is[j], is[j+1]
                    τ = system.t[a] - system.t[b]
                    Sᵢₒₙ += Eᵢₒₙₛ[j]*τ
                    Sₑₗ += volkov_phase(𝐩s[j], system.volkov, a, b)
                end
                S = Sᵢₒₙ + Sₑₗ
                aₚᵣₒₚ = prod(prefactors)*exp(-im*S)
            end

            verbosity > 10 && @show is prefactors

            ∂a = @timeit to "Prefactor" (ionization(system, diagram, 𝐩s[end], 𝐀, is[end]) *
                                         aₚᵣₒₚ *
                                         recombination(system, diagram, 𝐩s[1], 𝐀, iref))

            @timeit to "Interior interactions" begin
                # Loop over "interior" interactions
                for j = (order>1 && first(diagram)[2]==0 ? 3 : 2):order
                    ion,which = diagram.path[j-1]
                    α = ions[j-1]
                    β = ions[j]
                    verbosity > 20 && @show j, ion, which α,β
                    interaction = system.couplings[which][α,β]

                    𝐤ᵢ = 𝐩s[momenta[j-1]]
                    𝐩ᵢ = 𝐩s[momenta[j]]
                    𝐀ᵢ = 𝐀[is[j]]

                    ∂a *= @timeit to "Interaction" interaction(kinematic_momentum(𝐤ᵢ, 𝐀ᵢ), kinematic_momentum(𝐩ᵢ, 𝐀ᵢ), is[j+1])
                end
            end

            amplitude += ∂a
        end
    end

    @timeit to "Weighting" weight*amplitude
end

# * High-level interface

function photoelectron_spectrum(k::AbstractArray{T},
                                system::System, diagram::Diagram;
                                iref=length(system.t),
                                verbosity=1, kwargs...) where T
    verbosity > 0 && @info "Photoelectrum spectrum calculation" system diagram length(k)

    cT = complex(eltype(T))
    c = similar(k, cT)
    to = TimerOutput()
    p = Progress(length(k))
    @timeit to "k loop" begin
        threaded_range_loop(eachindex(k)) do i
            tok = TimerOutput()
            c[i] = integrate_diagram(cT, system, diagram, iref, k[i]; to=tok, verbosity=verbosity-1, kwargs...)
            ProgressMeter.next!(p)
            merge!(to, tok, tree_point=["k loop"])
        end
    end
    TimerOutputs.complement!(to)
    verbosity > 0 && print_timer(to)
    c
end

function photoelectron_spectrum(k, args...; kwargs...)
    system = System(args...; kwargs...)
    photoelectron_spectrum(k, system, Diagram(system); kwargs...)
end

"""
    induced_dipole(system, diagram[; kwargs...])

Compute the induced dipole moment as a function of time of the
[`System`](@ref), for a specific [`Diagram`](@ref).
"""
function induced_dipole(system::System, diagram::Diagram; verbosity = 1, kwargs...)
    verbosity > 0 && @info "Induced dipole calculation" system diagram

    t = system.t

    DT = eltype(system.𝐀)
    𝐝 = zeros(DT, length(t))

    memory = get(kwargs, :memory, typemax(Int))

    @showprogress for (i,t) in enumerate(t)
        𝐝̃ = integrate_diagram(DT, system, diagram, i; imin=i-memory, verbosity=verbosity-1, kwargs...)
        𝐝[i] = 2real(𝐝̃)
    end

    𝐝
end

function induced_dipole(args...; kwargs...)
    system = System(args...; kwargs...)
    diagram = Diagram([(1,0),(1,0)], system)
    induced_dipole(system, diagram; kwargs...)
end

# * Exports

export photoelectron_spectrum, induced_dipole
