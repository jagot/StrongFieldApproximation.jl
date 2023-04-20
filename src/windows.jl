flat_window(n=typemax(Int)) = Trues(n)

#=

Smooth turn-off of scalar functions, following Eqs. (18–21) of

- Becke, A. D. (1988). A Multicenter Numerical Integration Scheme for
  Polyatomic Molecules. The Journal of Chemical Physics, 88(4),
  2547–2553. http://dx.doi.org/10.1063/1.454033

=#

function becke_smoother(μ, k::Integer)
    p = μ -> 3μ/2 - μ^3/2
    f = μ
    for i = 1:k
        f = p(f)
    end
    (1-f)/2
end

function becke_smoother(x, a, b, k)
    if x < a
        one(x)
    elseif x > b
        zero(x)
    else
        becke_smoother(2*(x-a)/(b-a)-1, k)
    end
end

function smooth_window(::Type{T}, nflat, nsmooth, k) where T
    w = ones(T, nflat+nsmooth)
    for i = eachindex(w)
        w[i] = becke_smoother(i, nflat, nflat+nsmooth, k)
    end
    w
end

export flat_window, smooth_window
