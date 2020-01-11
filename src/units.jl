@derived_dimension ElectricField Unitful.𝐈^-1*Unitful.𝐋*Unitful.𝐌*Unitful.𝐓^-3
@derived_dimension Intensity Unitful.𝐌*Unitful.𝐓^-3
@derived_dimension InvIntensity Unitful.𝐌^-1*Unitful.𝐓^3

atomic_units(I::Intensity) = I/(3.5094452e16*u"W"/(u"cm"^2)) |> NoUnits
atomic_units(iI::InvIntensity) = iI*3.5094452e16*u"W"/(u"cm"^2) |> NoUnits
atomic_units(E::ElectricField) = E/(5.14220651e11*u"V"/u"m") |> NoUnits
atomic_units(q) = austrip(q)
