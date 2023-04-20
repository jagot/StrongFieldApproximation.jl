# * Diagrams

"""
    Diagram(path, couplings)

Goldstone diagram describing a strong-field process, i.e. after an
initial photoionization event, the ion and photoelectron are
propagated separately (at the chosen level of approximation, typically
Volkov waves are used for the photoelectrons). By interacting with the
external field, the ion state may change, and through electron
(re)scattering, both photoelectron and ion state may change.

The diagram is specified using a `path` of pairs (in antichronological
order), where each pair designates `(resultant ion channel,
interaction)`, where `interaction` is an index to one of the
`couplings`. The special value `interaction == 0` corresponds to the
initial ionization event, and should appear once at the end of `path`
(i.e. the first event); if it also appears as the first element of
`path` (i.e. the last event), it corresponds to recombination to the
initial state.

# Examples

```jldoctest
julia> couplings = (StrongFieldApproximation.DipoleCoupling,
                    StrongFieldApproximation.CoulombCoupling)
(StrongFieldApproximation.DipoleCoupling, StrongFieldApproximation.CoulombCoupling)

julia> Diagram([(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐤


julia> Diagram([(1,0),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ╲   ╱⇝
   |0⟩

julia> Diagram([(2,1),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
 ⇝┃   │
 2┃   │𝐤


julia> Diagram([(2,2),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ┠┈┈┈┤
 2┃   │𝐤


julia> Diagram([(1,2),(2,2),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ┠┈┈┈┤
 2┃   │𝐪
  ┠┈┈┈┤
 1┃   │𝐤


julia> Diagram([(3,1),(1,1),(2,2),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ┠┈┈┈┤
 2┃   │𝐪
 ⇝┃   │
 1┃   │
 ⇝┃   │
 3┃   │𝐤
```
"""
struct Diagram{Couplings}
    path::Vector{Tuple{Int,Int}}
    couplings::Couplings
    function Diagram(path, couplings::Couplings) where Couplings
        if !isempty(path)
            α,which = last(path)
            which == 0 ||
                throw(ArgumentError("Non-empty diagrams must have photoionization as the first interaction"))
        end
        ncoup = length(couplings)
        vcoup = 0:ncoup
        invalid_couplings = findall(∉(vcoup), last.(path))
        if !isempty(invalid_couplings)
            length(invalid_couplings) == 1 &&
                throw(ArgumentError("Coupling at vertex $(first(invalid_couplings)) not in valid range $(vcoup)"))
            throw(ArgumentError("Couplings at vertices $(invalid_couplings) not in valid range $(vcoup)"))
        end
        for i = 2:length(path)-1
            path[i][2] == 0 &&
                throw(ArgumentError("Only the first and last interaction in a diagram may have which = 0 (corresponding to ionization/recombination)"))
        end
        if length(path) > 1 && path[1][2] == 0
            path[1][1] == path[2][1] ||
                throw(ArgumentError("Must recombine from the last channel"))
        end
        new{Couplings}(path, couplings)
    end
end

Diagram(path::Tuple{Int,Int}, couplings) =
    Diagram([path], couplings)

"""
    Diagram(path::AbstractVector{Tuple{Int,Int}}, system)

Convenience constructor that sets up the [`Diagram`](@ref) given a
`path` and the couplings present in `system`. If `path` is empty, the
diagram corresponding to photoionization into the first ionization of
channel of `system` will automatically be generated.
"""
function Diagram(path::AbstractVector, system::System)
    nc = length(system.ionization_channels)
    if isempty(path)
        nc > 1 &&
            @warn "More than one ionization channel present in system, choosing the first"
        path = [(1,0)]
    end
    vc = 1:nc
    invalid_channels = findall(∉(vc), first.(path))
    if !isempty(invalid_channels)
        length(invalid_channels) == 1 &&
            throw(ArgumentError("Invalid channel at vertex $(first(invalid_channels))"))
        throw(ArgumentError("Invalid channels at vertices $(invalid_channels)"))
    end
    Diagram(path, [typeof(first(c)) for c in system.couplings])
end

Diagram(system::System) = Diagram([], system)

for f in [:length, :first, :firstindex, :lastindex, :isempty]
    @eval Base.$f(d::Diagram) = $f(d.path)
end
Base.getindex(d::Diagram, i) = Diagram(d.path[i], d.couplings)

# ** Pretty-printing

function draw_ionization(io, nd)
    println(io, lpad("|0⟩", nd+4))
    print(io, lpad("╱ ╲", nd+4))
    printstyled(io, "⇜", color=:light_red)
    println(io)
end

function draw_recombination(io, nd)
    print(io, lpad("╲ ╱", nd+4))
    printstyled(io, "⇝", color=:light_red)
    println(io)
    println(io, lpad("|0⟩", nd+4))
end

draw_exciton(io, ion, nd, electron="") =
    println(io, lpad(ion, nd)*"┃   │$(electron)")

draw_coulomb_interaction(io, nd) =
    println(io, repeat(" ",nd)*"┠┈┈┈┤")

function draw_dipole_interaction(io, nd)
    printstyled(io, lpad("⇝",nd), color=:light_red)
    println(io, "┃   │")
end

function Base.show(io::IO, ::MIME"text/plain", d::Diagram)
    println(io, "Goldstone Diagram:")
    if isempty(d)
        println(io, "|0⟩")
        return
    end
    nd = maximum(iw -> length(digits(iw[1])), d.path)+1
    electrons = ["𝐩","𝐪"]
    cur_electron = 1
    for (i,(ion,which)) in Iterators.reverse(enumerate(d.path))
        if which == 0
            if i == length(d.path)
                draw_ionization(io, nd)
            elseif i == 1
                draw_recombination(io, nd)
            end
        else
            c = d.couplings[which]
            if c <: CoulombCoupling
                draw_coulomb_interaction(io, nd)
            elseif c <: DipoleCoupling
                draw_dipole_interaction(io, nd)
            end
        end
        electron = if i == 1
            "𝐤"
        elseif !iszero(which) && d.couplings[which] <: DipoleCoupling
            ""
        else
            ei = mod1(cur_electron,length(electrons))
            e = electrons[ei]*repeat("′", fld1(cur_electron,length(electrons))-1)
            cur_electron += 1
            e
        end

        which == 0 && i == 1 && length(d) > 1 || draw_exciton(io, ion, nd, electron)
    end
end

# ** Diagram properties

function analyze_diagram(system, diagram)
    α,which = first(diagram)
    ions = Int[]
    # For photoelectron spectra, time 1 is the reference time,
    # typically at which the laser pulse has ended; for direct
    # photoelectron diagrams, time 2 is the time of ionization. For
    # dipoles, time 1 is the time of recombination.
    unique_momenta = Tuple{Int,Int}[(1,2)]
    momenta = [1]
    indeterminate_momenta = Int[]
    ld = length(diagram)
    order = ld

    if ld > 1
        if which == 0
            push!(indeterminate_momenta, 1)
            # Recombination does not increase the order of the diagram, since
            # it only amounts to projecting the wavefunction on ⟨0|𝐫.
            order -= 1
        end
    else
        push!(ions, α)
    end

    i = 2
    for (α,which) in diagram.path[(iszero(which) ? 2 : 1):end]
        push!(ions, α)
        a,b = unique_momenta[end]
        unique_momenta[end] = (a,i)
        if !(iszero(which) || canonical_momentum_conservation(system, which) == CanonicalMomentumConservation())
            push!(unique_momenta, (i,i+1))
            push!(indeterminate_momenta, length(unique_momenta))
        end
        if !iszero(which)
            push!(momenta, length(unique_momenta))
        end

        i += 1
    end

    return ions, unique_momenta, momenta, indeterminate_momenta, order
end

# 𝐤 is nothing, we want a dipole amplitude
function diagram_amplitude(::Type{Amp}, system::System, diagram::Diagram, iref, i, ::Nothing) where Amp
    amplitude = complex(zero(Amp))

    amplitude
end

# 𝐤 is determinate, we want a photoelectron spectrum
function diagram_amplitude(::Type{Amp}, system::System, diagram::Diagram, iref, i, 𝐤) where Amp
    amplitude = complex(zero(Amp))

    amplitude
end

# * Exports

export Diagram
