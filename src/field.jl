mutable struct AccumulatedIntegral{Tt,Ft,Func,Wt}
    t::Tt
    ∫f::Ft
    f::Func
    c::Vector{Tt} # Quadrature roots
    w::Vector{Wt} # Quadrature weights

    function AccumulatedIntegral(f::Func, ::Type{Ft}, ::Type{Tt}=Ft;
                                 k = 3,
                                 t::Tt=zero(Tt), init=zero(Ft)) where {Ft,Func,Tt}
        c,w = gausslegendre(k)
        new{Tt,Ft,Func,eltype(w)}(t, init, f, c, w)
    end
end

function (F::AccumulatedIntegral)(a,b)
    μ = (a+b)/2
    δ = (b-a)/2
    v = zero(F.∫f)
    for (c,w) in zip(F.c,F.w)
        v += w*F.f(δ*c + μ)
    end
    δ*v
end

(F::AccumulatedIntegral)(t) = F.∫f + F(F.t, t)

function set!(F::AccumulatedIntegral, t)
    F.∫f = F(t)
    F.t = t
end

function field(F, ndt; kwargs...)
    t = timeaxis(F, ndt)
    tau = atomic_units.(t)

    Ev = atomic_units.(F.(t))
    A = AccumulatedIntegral(t -> -atomic_units(F(24.2u"as"*t)), Float64;
                            t=first(tau), kwargs...)
    nt = length(t)
    Av = zeros(Float64, nt)
    for i = 1:nt
        Av[i] = A(tau[i])
        set!(A, tau[i])
    end

    tau,Ev,Av
end
