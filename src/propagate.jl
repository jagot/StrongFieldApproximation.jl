using SFA.trapz

#=
\[x(t)=
  i\int\limits_0^t d t'
  \int d^3\vec{p}\;
  d_x^*[\vec{p}-\vec{A}(t)]
  \exp[-i S(\vec{p},t,t')]
  E(t')
  d_x[\vec{p}-\vec{A}(t')]
  +\textrm{c.c.},\]
where
\[
  S(\vec{p},t,t')\equiv
  \int\limits_{t'}^t d t''
  \left\{
  \frac{[\vec{p}-\vec{A}(t'')]^2}{2}+I_p
  \right\}.
  \]

We approximate this using the saddle-point method:
\[x(t)\approx
  i\int\limits_0^\infty d\tau
  \left(\frac{\pi}{\varepsilon+i\tau/2}\right)^{3/2}
  d_x^*[p_\textrm{st}(t,\tau)-A_x(t)]
  \exp[-i S_\textrm{st}(t,\tau)]
  E(t-\tau)
  d_x[p_\textrm{st}(t,\tau)-A_x(t-\tau)]
  +\textrm{c.c.},\]
with \(\tau\equiv t-t'\), \(S_{\textrm{st}}(t,\tau)\equiv S(p_{\textrm{st}},t,t-\tau)\) and
\[p_{\textrm{st}}(t,\tau)=\frac{x(t)-x(t-\tau)}{\tau} = \frac{1}{\tau}\int\limits_{t-\tau}^t dt' A(t').\]
=#

function calc_i(i, offset, jmin, δt, Av, Ev, Ip, d)
    τ = reverse(1:i+offset)*δt
    nτ = length(τ)

    At = Av[i+offset]

    pτ = 0
    dx = 0im

    for j = nτ:-1:max(jmin,1)
        ζ = (π/(1e-10 + im*(τ[j]-δt)/2))^(3/2)

        AtP = Av[j]
        EtP = Ev[j]
        pτ += AtP
        p = pτ*δt/τ[j] # Saddle-point momentum

        pion = p - AtP
        prec = p - At

        dion = d(pion)
        drec = d(prec)

        S = trapzv(0.5(p .- Av[j:i+offset-1]).^2 .+ Ip, δt) # Saddle-point action

        dx += ζ*conj(drec)*exp(-im*S)*EtP*dion
    end
    -imag(dx)
end

function propagate(t::Vector,
                   A::Vector,E::Vector,
                   Ip::Real,d::Function,
                   T::Real,ndt::Integer,
                   jmin::Function, offset::Int)
    nt = length(t)
    δt = T/ndt

    x = Vector{Float64}(undef, nt)
    Threads.@threads for i = 1:nt
        x[i] = calc_i(i, offset, jmin(i), δt, A, E, Ip, d)
    end

    x *= 2δt
    x, t
end

function propagate(t::Vector,
                   A::Vector,E::Vector,
                   Ip::Real,d::Function,
                   T::Real,ndt::Integer,
                   tmin::Real = -0.65)
    imin = round(Int, tmin*ndt)
    jmin = if tmin < 0.0
        i -> i + imin
    else
        i -> imin
    end
    propagate(t, A, E, Ip, d,
              T, ndt, jmin, 0)
end

function propagate(F::ElectricFields.LinearField,
                   Ip::Real,d::Function,
                   ndt::Integer, args...;
                   kwargs...)
    tau,E,A = field(F, ndt; kwargs...)
    propagate(tau, A, E, Ip, d,
              atomic_units(period(F)), ndt, args...)
end

export propagate
