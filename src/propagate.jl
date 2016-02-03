using ProgressMeter
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

# tlims set the limits of the integration: t will assume values in the
# interval [tlims[1],tlims[2]] with ndt steps per cycle. t' will
# assume values in the interval [tlims[3],t], ∀ t, tlims[3] ⩾ 0, else
# [t - tlims[3],t].
function propagate(A,E,Ip,d,tlims,T,ndt)
    δt = T/ndt

    t = (tlims[1]*ndt+1:tlims[2]*ndt)*δt
    offset = (tlims[1]-tlims[3])*ndt-1
    nt = length(t)
    t2 = (2:nt+offset+1)*δt
    tlims = collect(tlims)*T

    Av = A(t2)
    Ev = E(t2)

    x = zeros(nt)

    @showprogress for i = nt:-1:1
        τ = if tlims[3] < 0
            [1]*δt # Not implemented yet
        else
            reverse(1:i+offset)*δt
        end
        nτ = length(τ)

        At = Av[i+offset]

        pτ = 0
        dx = 0im

        for j = nτ:-1:1
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
        x[i] = -imag(dx)
    end
    x *= 2δt

    x,t
end

export propagate
