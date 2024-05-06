struct ComplexPlane{T,R<:AbstractRange{T}} <: AbstractMatrix{Complex{T}}
    re::R
    im::R
end

Base.length(p::ComplexPlane) = length(p.re)*length(p.im)
Base.size(p::ComplexPlane) = (length(p.im), length(p.re))
Base.eltype(p::ComplexPlane) = complex(eltype(p.re))
Base.getindex(p::ComplexPlane, i::CartesianIndex) = p.re[i[2]] + im*p.im[i[1]]

function primitive_integral(f, a, b, c, w)
    μ = (a+b)/2
    δ = (b-a)/2
    v = zero(f(a))
    for (c,w) in zip(c,w)
        v += w*f(δ*c + μ)
    end
    δ*v
end

interval(x::AbstractRange, r::UnitRange) = x[r[1]]..x[r[end]]

function split_range(r::UnitRange)
    n = length(r)
    n < 1 && error("Too small range")
    n == 1 && return (r, r[1]:(r[1]-1))

    s = n ÷ 2
    r[1]:r[s], r[s+1]:r[end]
end

function find_closest(x::AbstractRange, x₀::Real)
    cur = firstindex(x):lastindex(x)
    i = 0
    while length(cur) > 1 && i < ceil(Int, log2(length(x)))
        left,right = split_range(cur)
        if x₀ ∈ interval(x, left)
            cur = left
            continue
        elseif x₀ ∈ interval(x, right)
            cur = right
            continue
        end
        x₀ < x[left[1]] && return left[1]
        x₀ > x[right[end]] && return right[end]
        return norm(x[left[end]]-x₀) < norm(x[right[1]]-x₀) ? left[end] : right[1]
        i += 1
    end
    cur[1]
end

function find_closest(x::ComplexPlane, x₀::Number)
    j = find_closest(x.re, real(x₀))
    i = find_closest(x.im, imag(x₀))
    CartesianIndex((i,j))
end

struct AutomaticIntegration{Func, Grid, Xt, GridIntegrals,
                            Roots, Weights}
    f::Func
    x::Grid
    x₀::Xt
    ∫f::GridIntegrals
    c::Roots
    w::Weights
end

function integrate!(ai::AutomaticIntegration{<:Any,<:AbstractRange}, i)
    a,b = firstindex(ai.x),lastindex(ai.x)
    for j = i+1:b
        ai.∫f[j] = ai(ai.x[j],i=j-1)
    end
    for j = i-1:-1:a
        ai.∫f[j] = ai(ai.x[j],i=j+1)
    end
end

function integrate!(ai::AutomaticIntegration{<:Any,<:ComplexPlane}, i)
    s = size(ai.x)
    ci = CartesianIndices(s)

    north(I::CartesianIndex) = CartesianIndex(I[1]-1,I[2])
    south(I::CartesianIndex) = CartesianIndex(I[1]+1,I[2])
    west(I::CartesianIndex) = CartesianIndex(I[1],I[2]-1)
    east(I::CartesianIndex) = CartesianIndex(I[1],I[2]+1)

    # First, we integrate parallel to the real axis, starting from
    # index i, then we integrate parallel to the imaginary axis,
    # starting from the real line-parallel we just created.
    for (r,dir) in [(reverse(ci[i[1],1:i[2]-1]), east),
                    (ci[i[1],i[2]+1:end], west),
                    (ci[i[1]+1:end,:], north),
                    (ci[i[1]-1:-1:1,:], south)]
        for I in r
            # @show I, dir(I) ai.∫f[I], ai.∫f[dir(I)]
            ai.∫f[I] = ai(ai.x[I],i=dir(I))
        end
    end
end

function AutomaticIntegration(f, x, x₀; k=3, init=zero(f(x₀)))
    Xt = eltype(x)
    s = size(x)
    Ft = typeof(f(Xt(x₀)))
    ∫f = zeros(Ft, s)

    c,w = gausslegendre(k)
    ai = AutomaticIntegration(f, x, x₀, ∫f, c, w)

    # We start by integrating from the point on the grid closest to x₀
    # to x₀, which is minus the integral value desired, since the
    # integral originates at x₀.
    i = find_closest(x, x₀)
    ∫f[i] = init - ai(x₀, i=i)
    integrate!(ai, i)

    ai
end

(ai::AutomaticIntegration)(x; i=find_closest(ai.x, x)) =
    ai.∫f[i] + primitive_integral(ai.f, ai.x[i], x, ai.c, ai.w)

(ai::AutomaticIntegration)(a, b) = ai(b) - ai(a)
