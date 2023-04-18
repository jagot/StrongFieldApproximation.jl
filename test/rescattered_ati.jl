# This rather unphysical combination of electric field and vector
# potential (disregarding short-pulse effects, i.e. the chain rule) is
# taken from Eqs. (21–22) of
#
# - D B Milo\vsevi\'c, Paulus, G. G., Bauer, D., & Becker,
#   W. (2006). Above-Threshold Ionization By Few-Cycle Pulses. Journal
#   of Physics B: Atomic, Molecular and Optical Physics, 39(14),
#   203–262. http://dx.doi.org/10.1088/0953-4075/39/14/r01
#
# to be able to reproduce the results shown in their Fig. 14.
function milosevic_pulse(F₀, ω, cycles, ndt, ϕ)
    T = 2π/ω
    steps = ceil(Int, cycles*ndt)
    t = range(0, stop=cycles*T, length=steps)
    Fv = zeros(steps)
    Av = zeros(steps)

    ωₚ = ω/cycles
    ω₀ = ω
    ω₁ = ω + ωₚ
    ω₂ = ω - ωₚ

    ϕ₁ = ϕ + π/2

    ℰ₀ = F₀/2
    ℰ₁ = -ℰ₀/2
    ℰ₂ = -ℰ₀/2

    𝒜₀ = ℰ₀/ω₀
    𝒜₁ = ℰ₁/ω₁
    𝒜₂ = ℰ₂/ω₂

    for (i,t) in enumerate(t)
        φ₀₁ = ω₀*t + ϕ₁
        φ₁₁ = ω₁*t + ϕ₁
        φ₂₁ = ω₂*t + ϕ₁

        Fv[i] = ℰ₀*sin(φ₀₁) + ℰ₁*sin(φ₁₁) + ℰ₂*sin(φ₂₁)
        Av[i] = 𝒜₀*cos(φ₀₁) + 𝒜₁*cos(φ₁₁) + 𝒜₂*cos(φ₂₁)
    end

    t, Fv, Av
end

@testset "Rescattered ATI" begin
    kmax = 3.0
    nk = 200
    nθ = 2
    spacing=[:momentum,:energy][2]
    𝐤,kmag,θ = momentum_grid(0.1, kmax, nk, nθ,
                            spacing=spacing)

    @field(F) do
        λ = 800.0u"nm"
        I₀ = 1e14u"W/cm^2"
        cycles = 4.0
        ϕ = π
        env = :cos²
    end

    ndt = ceil(Int, 3*(kmax^2/2)/photon_energy(F))

    t, Fv, Av = milosevic_pulse(amplitude(F), photon_energy(F), F.env.cycles, ndt, 0)

    Iₚ = 14u"eV" # "Krypton"

    # Elastic scattering off a Yukawa potential
    cc = StrongFieldApproximation.CoulombCoupling((𝐤,𝐩) -> yukawa_fourier(𝐩-𝐤, 1, 0, 1))
    couplings=[[cc;;]]

    ar = (Fv, Av, t)

    channel = IonizationChannel(Iₚ, ar...)
    system = StrongFieldApproximation.System([channel], nothing, couplings, ar...)

    direct_diagram = Diagram(system)
    rescattering_diagram = Diagram([(1,1),(1,0)], system)

    cdirect_ref = readdlm(reffile("krypton-ref-direct.csv"), ',', ComplexF64)
    crescattered_ref = readdlm(reffile("krypton-ref-rescattered.csv"), ',', ComplexF64)

    cdirect = photoelectron_spectrum(𝐤, system, direct_diagram, verbosity=1)
    test_approx_eq(cdirect, cdirect_ref)
    cdirect2 = multi_channel_photoelectron_spectrum(𝐤, system, direct_diagram, verbosity=1)
    @test size(cdirect2, 1) == 1
    test_approx_eq(cdirect2[1,:,:], cdirect_ref)

    crescattered = photoelectron_spectrum(𝐤, system, rescattering_diagram, verbosity=1)
    test_approx_eq(crescattered, crescattered_ref)
    crescattered2 = multi_channel_photoelectron_spectrum(𝐤, system, rescattering_diagram, verbosity=1)
    @test size(crescattered2, 1) == 1
    test_approx_eq(crescattered2[1,:,:], crescattered_ref)
end
