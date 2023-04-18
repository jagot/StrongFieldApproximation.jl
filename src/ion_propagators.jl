abstract type IonPropagator end

function interaction(ions::IonPropagator, coupling, α, k, β, p, i)
    s = zero(eltype(ions))
    for (j′,q′) in non_zero_ion_mapping(ions, α, i)
        for (j,q) in non_zero_ion_mapping(ions, β, i)
            s += conj(q′)*q*coupling[j′,j](k, p, i)
        end
    end
    s
end

function apply_interaction!(v, coupling, u, k, p, ti)
    nc = length(u)
    v .= false
    for α = 1:nc
        for β = 1:nc
            v[α] += coupling[α,β](k, p, ti)*u[β]
        end
    end
    v
end

function apply_interaction!(v, ions::IonPropagator, coupling, u, k, p, i)
    Q = ion_mapping(ions, i)
    mul!(v, Q, u)
    apply_interaction!(u, coupling, v, k, p, i)
    mul!(v, Q, u)
end

# * Laser-dressed ions

@doc raw"""
    LaserDressedIons

Pre-propagated laser-dressed basis for the ion states, i.e. the
time-dependent eigenbasis of the ionic Hamiltonian including the
external electric field. In the eigenbasis, the ionic propagator is
diagonal, and the entries are

```math
\matrixel{\eta}{\propU(a,b)}{\gamma} =
\exp[-\im E_\eta(T-a)]
\exp[+\im E_\eta(T-b)]
\delta_{\eta\gamma},
```

i.e. the reference time is ``T`` (the end of the pulse), which is
consistent with [`VolkovPhases`](@ref).
"""
struct LaserDressedIons{T,IonBasis} <: IonPropagator
    expimE::Matrix{T}
    ion_basis::IonBasis
end

function LaserDressedIons(E::Matrix{T}, Q, t) where T
    expimE = similar(E, complex(T))

    nt = length(t)
    for i = 1:nt
        μ = im*(t[end] - t[i])
        expimE[i,:] = exp.(μ.*E[i,:])
    end

    LaserDressedIons(expimE, Q)
end

function LaserDressedIons(ics, 𝐫ᵢₒₙ::SparseMatrixCSC, Fv::AbstractArray, t)
    to = TimerOutput()

    @timeit to "Allocations" begin
        m,n = size(𝐫ᵢₒₙ)
        @assert m == n == length(ics)

        nt = length(t)
        E₀ = [ic.E for ic in ics]
        H₀ = Diagonal(E₀)

        X = .!iszero.(I + (𝐫ᵢₒₙ .≠ 0))
        bs = find_blocks(X)

        @info "Diagonalizing $(m)×$(n) ionic Hamiltonian for $(nt) times" bs

        T = eltype(E₀)
        E = zeros(T, nt, m)
        Q = zeros(T, m, m, nt)
    end

    @timeit to "Block loop" for b in bs
        nb = length(b)
        if nb == 1
            j = b[1]
            E[:,j] .= E₀[j]
            Q[j,j,:] .= 1
            continue
        end

        Hsub = zeros(T, nb, nb)
        H₀sub = Diagonal(E₀[b])
        zsub = Matrix(𝐫ᵢₒₙ[b,b])
        @timeit to "Time loop" begin
            @showprogress for i = 1:nt
                Hsub .= H₀sub .+ Fv[i] .* zsub
                @assert Hsub ≈ Hsub'
                ee = eigen!(Hermitian(Hsub))
                for (ij,j) in enumerate(b)
                    E[i,j] = ee.values[ij]
                end
                Q[b,b,i] = ee.vectors
            end
        end
    end

    TimerOutputs.complement!(to)
    print_timer(to)

    LaserDressedIons(E, Q, t)
end

function LaserDressedIons(ics, ::Nothing, _::AbstractArray, t)
    nt = length(t)
    E₀ = [ic.E for ic in ics]
    nc = length(E₀)

    T = eltype(E₀)
    E = zeros(T, nt, nc)
    for j = 1:nc
        E[:,j] .= E₀[j]
    end

    LaserDressedIons(E, nothing, t)
end

LaserDressedIons(ics, 𝐫ᵢₒₙ, F::ElectricFields.AbstractField, t) =
    LaserDressedIons(ics, 𝐫ᵢₒₙ, field_amplitude(F, t), t)

Base.eltype(::LaserDressedIons{T}) where T = T

ion_mapping(ions::LaserDressedIons, i) =
    view(ions.ion_basis, :, :, i)

ion_mapping(ions::LaserDressedIons{<:Any,Nothing}, _) = I

ion_mapping(ions::LaserDressedIons, α, i) =
    view(ions.ion_basis, α, :, i)

ion_mapping(ions::LaserDressedIons{T,Nothing}, α, _) where T =
    vcat(zeros(T,α-1), one(T), zeros(T, size(ions.expimE,2)-α))

non_zero_ion_mapping(ions::LaserDressedIons, α, i) =
    filter(jq -> !iszero(last(jq)),
           collect(enumerate(ion_mapping(ions, α, i))))

non_zero_ion_mapping(ions::LaserDressedIons{T,Nothing}, α, _) where T =
    [(α, one(T))]

ion_propagation(ions::LaserDressedIons, j, a, b) =
    ions.expimE[a,j] .* conj(ions.expimE[b,j])

# * Field-free ions

struct FieldFreeIons{T,Times} <: IonPropagator
    E₀::Vector{T}
    t::Times
end

FieldFreeIons(ics, _, _, t) =
    FieldFreeIons([ic.E for ic in ics], t)

Base.eltype(::FieldFreeIons{T}) where T = T

ion_mapping(::FieldFreeIons, _) = I

ion_mapping(::FieldFreeIons{T}, α, _) where T =
    vcat(zeros(T,α-1), one(T), zeros(T, length(ions.E₀)-α))

non_zero_ion_mapping(::FieldFreeIons{T}, α, _) where T =
    [(α, one(T))]

ion_propagation(ions::FieldFreeIons, j, a, b) =
    exp.(-im*(ions.t[a]-ions.t[b])*ions.E₀[j])
