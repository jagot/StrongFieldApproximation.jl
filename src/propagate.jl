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

    for j = nτ:-1:jmin
        ζ = (π/(1e-10 + im*(τ[j]-δt)/2))^(3/2)

        AtP = Av[j]
        EtP = Ev[j]
        pτ += AtP
        p = pτ*δt/τ[j] # Saddle-point momentum

        pion = p - AtP
        prec = p - At

        dion = d(pion)
        drec = d(prec)

        S = trapzv(0.5(p - Av[j:i+offset-1]).^2 + Ip, δt) # Saddle-point action

        dx += ζ*conj(drec)*exp(-im*S)*EtP*dion
    end
    -imag(dx)
end

# tlims set the limits of the integration: t will assume values in the
# interval [tlims[1],tlims[2]] with ndt steps per cycle. t' will
# assume values in the interval [tlims[3],t], ∀ t, tlims[3] ⩾ 0, else
# [t - tlims[3],t].
function propagate(A::Function,E::Function,Ip::Real,d::Function,tlims,T::Real,ndt::Integer)
    δt = T/ndt

    t = (tlims[1]*ndt+1:tlims[2]*ndt)*δt
    if tlims[3]>=0
        offset = (tlims[1]-tlims[3])*ndt-1
        jmin = 1
    else
        offset = tlims[1]*ndt - 1
        jmin = round(Int, offset + tlims[3]*ndt)
    end
    nt = length(t)
    tlims = collect(tlims)*T

    t2 = SharedArray(Float64, nt+offset)
    Av = SharedArray(Float64, nt+offset)
    Ev = SharedArray(Float64, nt+offset)
    @parallel for i = 1:nt+offset
        tt = (i+1)*δt
        t2[i] = tt
        Av[i] = A(tt)
        Ev[i] = E(tt)
    end

    x = sdata(SharedArray(Float64, nt,
                          init = S -> (S[Base.localindexes(S)] =
                                      map(i -> calc_i(i, offset, jmin + i, δt, Av, Ev, Ip, d),
                                          Base.localindexes(S)))))
    x *= 2δt

    x,t
end

function propagate(t::AbstractVector,
                   A::Vector,E::Vector,
                   Ip::Real,d::Function,
                   T::Real,ndt::Integer,
                   tmin::Real = -0.65)
    nt = length(t)

    Av = SharedArray(Float64, nt)
    Ev = SharedArray(Float64, nt)
    @parallel for i = 1:nt
        Av[i] = A[i]
        Ev[i] = E[i]
    end

    δt = T/ndt
    jmin = round(Int, tmin*ndt)

    x = sdata(SharedArray(Float64, nt,
                           init = S -> (S[Base.localindexes(S)] =
                                       map(i -> SFA.calc_i(i, 0, max(jmin+i, 1), δt, Av, Ev, Ip, d),
                                           Base.localindexes(S)))))

    x *= 2δt
    x
end

export propagate
