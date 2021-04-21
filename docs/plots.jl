using StrongFieldApproximation
using ElectricFields
using Unitful
using UnitfulAtomic
using FFTW

using PyPlot
using Jagot.plotting
plot_style("ggplot")

using Jagot

function savedocfig(name,dir="figures")
    fig = gcf()
    filename = joinpath(@__DIR__, "src", dir, "$(name).svg")
    savefig(filename,
            transparent=false,
            facecolor=fig.get_facecolor())
    # close(fig)
    if isfile(filename)
        println("Saved $(name) to $(filename)")
    else
        @warn "Saving $(name) to $(filename) failed"
    end
end

function savedfigure(fun::Function, name, args...; kwargs...)
    cfigure(name, args...; kwargs...) do
        fun()
    end
    savedocfig(name)
end

spectrum(d::AbstractVector{<:Number}) = fftshift(fft(d))

function spectrum(d::AbstractVector{<:AbstractVector{T}}) where {T}
    D = Matrix{complex(T)}(undef, length(d), 3)
    for (i,dd) in enumerate(d)
        D[i,:] .= dd
    end
    fft!(D, 1)
    fftshift(D, 1)
end

function hhg_example()
    @field(F) do # Linear polarization
        λ = 800.0u"nm"
        I₀ = 3e14u"W/cm^2"
        τ = 6.2u"fs"
        toff = 10.0u"fs"
        tmax = 13.0u"fs"
        env = :trunc_gauss
    end;

    @field(F2) do # Slightly elliptical polarization
        λ = 800.0u"nm"
        I₀ = 3e14u"W/cm^2"
        τ = 6.2u"fs"
        toff = 10.0u"fs"
        tmax = 13.0u"fs"
        env = :trunc_gauss
        ξ = 0.2
    end;

    ndt = 300 # Steps per cycle
    Iₚ = 0.5 # Hydrogen

    # d will be a vector of scalars
    system,diagram,d = induced_dipole(Iₚ, F, ndt, memory=floor(Int, 0.65ndt));

    # d2 will be a vector of 3d vectors
    system2,diagram2,d2 = induced_dipole(Iₚ, F2, ndt, memory=floor(Int, 0.65ndt));

    t = timeaxis(F, ndt)
    tplot = 24.2e-3t
    q = fftshift(fftfreq(length(t), ndt))
    qsel = ind(q,0):ind(q,100)
    D = spectrum(d)
    D2 = spectrum(d2)
    cutoff = 3.17austrip(ponderomotive_potential(F)) + Iₚ


    savedfigure("hhg_example", figsize=(8,8)) do
        csubplot(311) do
            plot(tplot, field_amplitude(F, t))
            plot(tplot, field_amplitude(F2, t), "--")
            legend(["1d", "3d x", "3d y", "3d z"])
            axes_labels_opposite(:x)
            xlabel(L"$t$ [fs]")
            ylabel(L"$F(t)$ [au]")
        end
        csubplot(312,nox=true) do
            plot(tplot, d)
            plot(tplot, d2, "--")
            legend(["1d", "3d x", "3d y", "3d z"])
            ylabel(L"\langle\mathbf{r}\rangle(t)")
        end
        csubplot(313) do
            semilogy(q[qsel], abs2.(D[qsel]))
            semilogy(q[qsel], abs2.(D2[qsel,:]), "--")
            legend(["1d", "3d x", "3d y", "3d z"])
            axvline(cutoff/photon_energy(F), linestyle=":", color="black")
            xlabel(L"Harmonic order of 800 nm [$q$]")
            ylabel(L"|\mathbf{r}(q)|^2")
        end
    end
end

function streaking_example()
    @field(pump) do
        ħω = 80u"eV"
        I₀ = 1e8u"W/cm^2"
        τ = 200u"as"
        toff = 400u"as"
        tmax = 500u"as"
        env = :trunc_gauss
    end;

    @field(probe) do
        λ = 800.0u"nm"
        I₀ = 1e12u"W/cm^2"
        τ = 2.66u"fs"
        toff = 4.0u"fs"
        tmax = 5.0u"fs"
        env = :trunc_gauss
    end;

    F = pump + delay(probe, 1.0u"fs");

    ndt = 100

    nk = 100
    nθ = 2 # Forward and backward spectra only
    k,kmag,θ = momentum_grid(50u"eV", 80u"eV", nk, nθ,
                             spacing=:energy);

    Iₚ = 14u"eV" # "Krypton"
    c = photoelectron_spectrum(k, Iₚ, F, ndt);

    x = 27.211kmag.^2/2

    savedfigure("streaking_single_delay") do
        plot(x, abs2.(c))
        legend(["Forward", "Backward"])
        xlabel(L"$W_k$ [eV]")
    end

    delays = range(-4u"fs", stop=4u"fs", length=100)

    forward = zeros(ComplexF64, nk, length(delays))
    backward = zeros(ComplexF64, nk, length(delays))

    for (j,d) in enumerate(delays)
        println(j,"/",length(delays))
        cc = photoelectron_spectrum(k, Iₚ, pump + delay(probe, d), ndt, verbosity=0)
        forward[:,j] = cc[:,1]
        backward[:,j] = cc[:,2]
    end

    savedfigure("streaking_spectrogram", figsize=(8,10)) do
        csubplot(311) do
            plot(ustrip(delays), vector_potential(probe, austrip.(delays)))
            axes_labels_opposite(:x)
            xlabel(L"$t$ [fs]")
            ylabel(L"$A_{\mathrm{IR}}(t)$ [au]")
            margins(0,0.1)
        end
        csubplot(312, nox=true) do
            plot_map(ustrip(delays), x, abs2.(forward))
            ylabel(L"$W_k$ forward [eV]")
        end
        csubplot(313) do
            plot_map(ustrip(delays), -x, abs2.(backward))
            xlabel("Pump–probe delay [fs]")
            ylabel(L"$W_k$ backward [eV]")
        end
    end

end

function ppt_example()
    I = 10 .^ range(10, stop=15, length=100) * u"W/cm^2";
    Iau = ElectricFields.Iaustrip.(I);
    λ = [2000u"nm", 800u"nm", 500u"nm"]
    ω = austrip.(2π*1u"c" ./ λ)

    Iₚ = 0.5 # Hydrogen

    rates = StrongFieldApproximation.IonizationRates.PPT.(Iₚ, Iau, ω', 0, 0);

    savedfigure("ppt_example") do
        loglog(ustrip.(I), rates)
        xlabel(L"$I$ [W/cm$^2$]")
        ylabel(L"Ionization rates [s$^{-1}$]")
        legend(string.(λ))
    end
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
mkpath("docs/src/figures")
@echo hhg_example()
@echo streaking_example()
@echo ppt_example()
