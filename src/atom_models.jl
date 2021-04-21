#=
Gaussian model:
\[\vec{d}(\vec{p}) = -\mathrm{i}\left(\frac{1}{\pi\alpha}\right)^{3/4}
\frac{\vec{p}}{\alpha}e^{-p^2/2\alpha}\]
=#
d_gaussian(p,α) = -im*(1/(π*α))^(3/4) * p/α*exp(p^2/(2α))

#=
Hydrogenlike atom:
\[\vec{d}(\vec{p}) = -\mathrm{i}\left(\frac{2^{7/2}\alpha^{5/4}}{\pi}\right)
\frac{\vec{p}}{(p^2+\alpha)^3}\]
=#
d_hyd(p,α) = -im*(2^(7/2)*α^(5/4)/π)*p/((norm(p)^2+α)^3)
hyd_α(Ip) = 2Ip

hydrogen_like(Ip) = Ip, p -> d_hyd(p, hyd_α(Ip))
hydrogen = hydrogen_like(0.5)

function yukawa_fourier(q, a, b, λ)
    # Eq. (47) of
    #
    # - D B Milo\vsevi\'c, Paulus, G. G., Bauer, D., & Becker,
    #   W. (2006). Above-Threshold Ionization By Few-Cycle Pulses. Journal
    #   of Physics B: Atomic, Molecular and Optical Physics, 39(14),
    #   203–262. http://dx.doi.org/10.1088/0953-4075/39/14/r01
    C = norm(q)^2 + λ^2
    -(2b*λ + a*C)/(2π^2*C^2)
end

coulomb_fourier(q, Z=1) = yukawa_fourier(q, Z, 0, 0)

export d_gaussian, d_hyd, hyd_α, hydrogen_like, hydrogen,
    yukawa_fourier, coulomb_fourier
