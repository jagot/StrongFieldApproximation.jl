This implements the strong-field ionization formula given by PPT, i.e.

- Perelomov, A., Popov, V., & Terent'ev, M. (1966). Ionization of
  Atoms in an Alternating Electric Field. Soviet Physics—Journal of
  Experimental and Theoretical Physics, [23(5),
  924–934](http://www.jetp.ac.ru/cgi-bin/e/index/e/23/5/p924?a=list),

with some sprinkles of ADK:

- Ammosov, M. V., Delone, N. B., & Krainov, V. P. (1986). Tunnel
  Ionization of Complex Atoms and of Atomic Ions in Alternating
  Electromagnetic Field. Soviet Physics—Journal of Experimental and
  Theoretical Physics, [64(6),
  1191–1194](http://www.jetp.ac.ru/cgi-bin/e/index/e/64/6/p1191?a=list).

```@docs
StrongFieldApproximation.IonizationRates.PPT
```

```jldoctest
julia> using StrongFieldApproximation, ElectricFields, Unitful, UnitfulAtomic

julia> I = 10 .^ range(10, stop=15, length=100) * u"W/cm^2";

julia> Iau = ElectricFields.Iaustrip.(I);

julia> λ = [400u"nm", 800u"nm", 2000u"nm"]
3-element Vector{Quantity{Int64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
  400 nm
  800 nm
 2000 nm

julia> ω = austrip.(2π*1u"c" ./ λ)
3-element Vector{Float64}:
 0.11390838132237764
 0.05695419066118882
 0.022781676264475526

julia> Iₚ = 0.5 # Hydrogen
0.5

julia> rates = StrongFieldApproximation.IonizationRates.PPT.(Iₚ, Iau, ω', 0, 0);
```

![PPT example](figures/ppt_example.svg)

# Internals reference

```@autodocs
Modules = [StrongFieldApproximation.IonizationRates]
Filter = !isequal(StrongFieldApproximation.IonizationRates.PPT)
```
