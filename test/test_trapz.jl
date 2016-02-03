using PyCall
pygui(:qt)
using PyPlot

using SFA.trapz

using Calculus

FFs = [x -> x^4 - 4x^2, x -> x, x -> sin(2π*x).^3, x -> sin(x), x -> cos(x)]

for j = 1:length(FFs)
    FF = FFs[j]
    f = FF'
    xmax = 5
    x = linspace(0,xmax,xmax*100)
    t = Trapz(f,x)
    dx = (x[2]-x[1])

    F = zeros(x)
    @time for i = 1:length(x)
        F[i] = next!(t)
    end

    figure(j)
    clf()
    plot(x,map(f, x), "-", label=L"$f(x)=F'(x)$")
    plot(x,F, "-", label=L"$\int_a^b f(x)dx$, trapezoidal rule")
    plot(x,map(FF, x) - FF(x[1]), "--", label=L"$F(x)-F(a)$, analytical")
    xlabel(L"$x$")
    legend(framealpha=0.75)
    margins(0,0.1)
    tight_layout()
end

show()
