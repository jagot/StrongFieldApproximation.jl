#=
Gaussian model:
\[\vec{d}(\vec{p}) = i\left(\frac{1}{\pi\alpha}\right)^{3/4}
\frac{\vec{p}}{\alpha}e^{-p^2/2\alpha}\]
=#
d_gaussian(p,α) = im*(1/(π*α))^(3/4) * p/α.*exp(p.^2/(2α))

#=
Hydrogenlike atom:
\[\vec{d}(\vec{p}) = i\left(\frac{2^{7/2}\alpha^{5/4}}{\pi}\right)
\frac{\vec{p}}{(p^2+\alpha)^3}\]
=#
d_hyd(p,α) = im*(2^(7/2)*α^(5/4)/π)*p./((p.^2+α).^3)
hyd_α(Ip) = 2Ip

hydrogen_like(Ip) = Ip, p -> d_hyd(p, hyd_α(Ip))
hydrogen = hydrogen_like(0.5)

export d_gaussian, d_hyd, hyd_α, hydrogen_like, hydrogen
