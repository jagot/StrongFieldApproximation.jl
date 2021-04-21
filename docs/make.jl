using StrongFieldApproximation
using Documenter

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(;
    modules=[StrongFieldApproximation],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    repo="https://github.com/jagot/StrongFieldApproximation.jl/blob/{commit}{path}#{line}",
    sitename="StrongFieldApproximation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://www.tipota.org/StrongFieldApproximation.jl",
        assets=String[],
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict(
                    :defd => "≝",
                    :abs => ["\\left|#1\\right|", 1],
                    :arcsinh => "\\operatorname{arcsinh}",
                    :ket => ["|#1\\rangle", 1],
                    :bra => ["\\langle#1|", 1],
                    :braket => ["\\langle#1|#2\\rangle", 2],
                    :ketbra => ["|#1\\rangle\\!\\langle#2|", 2],
                    :matrixel => ["\\langle#1|#2|#3\\rangle", 3],
                    :vec => ["\\mathbf{#1}", 1],
                    :mat => ["\\mathsf{#1}", 1],
                    :conj => ["#1^*", 1],
                    :im => "\\mathrm{i}",
                    :operator => ["\\mathfrak{#1}", 1],
                    :Hamiltonian => "\\operator{H}",
                    :hamiltonian => "\\operator{h}",
                    :Lagrangian => "\\operator{L}",
                    :fock => "\\operator{f}",
                    :lagrange => ["\\epsilon_{#1}", 1],
                    :vary => ["\\delta_{#1}", 1],
                    :onebody => ["(#1|#2)", 2],
                    :twobody => ["[#1|#2]", 2],
                    :twobodydx => ["[#1||#2]", 2],
                    :direct => ["{\\operator{J}_{#1}}", 1],
                    :exchange => ["{\\operator{K}_{#1}}", 1],
                    :diff => ["\\mathrm{d}#1\\,", 1],
                    :saddle => ["#1^{(\\textrm{st})}", 1],
                    :ce => "\\mathrm{e}"
                )
            )
        ))
    ),
    pages=[
        "Home" => "index.md",
        "Strong-field ionization" => "strong_field_ionization.md",
        "Reference" => "reference.md",
    ],
    doctest = false
)

deploydocs(;
    repo="github.com/jagot/StrongFieldApproximation.jl",
)
