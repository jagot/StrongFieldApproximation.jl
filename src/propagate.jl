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

function calc_i(i, j, δt, Av, Ev, Ip, d, window)
    # By reversing the τ vector, we make sure that the smallest
    # excursion times are paired with the latest ionization events (t'
    # close to t).
    τ = reverse(1+j-j[1])*δt
    nτ = length(τ)

    At = Av[j[end]]

    pτ = 0
    dx = 0im

    # The loop goes from the latest ionization event before t until
    # the earliest.
    for jj = reverse(1:length(j))
        ζ = (π/(1e-10 + im*(τ[jj]-δt)/2))^(3/2)

        jjj = j[jj]

        AtP = Av[jjj]
        EtP = Ev[jjj]
        pτ += AtP
        p = pτ*δt/τ[jj] # Saddle-point momentum

        pion = p - AtP
        prec = p - At

        dion = d(pion)
        drec = d(prec)

        S = trapzv(0.5(p - Av[jjj:j[end]]).^2 + Ip, δt) # Saddle-point action

        dx += ζ*conj(drec)*exp(-im*S)*EtP*dion*window(i,jjj)[1]
    end
    -imag(dx)
end

function limits(tmin,tmax,τmax,ndt,windowed=false)
    i = 1:(tmax-tmin)*ndt # Index vector of output
    if τmax == 0
        τmax2 = 0
        j = i -> 1:i
        window = (i,jj) -> ones(Float64,length(jj))
    elseif τmax < 0
        τmax2 = τmax - (windowed ? 0.35 : 0)
        jmin = round(Int, abs(τmax2)*ndt)
        j = i -> i + (0:jmin-1)
        window = if !windowed
            (i,jj) -> ones(Float64,length(jj))
        else
            # Optional, soft window function used when selecting excursion
            # times < abs(τmax).
            w = j -> (j .< 0) + (1.0 - exp(-((j/ndt-1).^6/0.001)).^4) .* ((0.<=j) .* (j.<ndt))
            (i,jj) -> w(ndt+i-jj)
        end
    else
        error("Invalid τmax, $(τmax)")
    end
    # Real time vector, x will be calculated at these points.
    t = (tmin*ndt+1:tmax*ndt)/ndt
    # Time vector, including excursions before tmin, in the case of
    # τmax ≠ 0. Fields are to be evaluated at these points.
    tt = ((tmin+τmax2)*ndt+1:tmax*ndt)/ndt


    i,j,t,tt,window
end

# tlims set the limits of the integration: t will assume values in the
# interval [tlims[1],tlims[2]] with ndt steps per cycle. t' will
# assume values in the interval [t - |tlims[3]|,t].
function propagate(A::Function,E::Function,
                   Ip::Real,d::Function,
                   tlims,T::Real,ndt::Integer;
                   windowed = false)
    δt = T/ndt
    i,j,t,tt,window = limits(tlims...,ndt,windowed)
    nt = length(t)
    ntt = length(tt)

    Av = Array(Float64, ntt)
    Ev = Array(Float64, ntt)
    for (i,ttt) in enumerate(tt)
        Av[i] = A(ttt*ndt*δt)
        Ev[i] = E(ttt*ndt*δt)
    end

    x = Array(Float64, nt)
    for i in eachindex(t)
        x[i] = calc_i(i, j(i), δt, Av, Ev, Ip, d, window)
    end

    x*2δt,t*ndt*δt
end

export propagate
