using StrongFieldApproximation
using DSP
using PyPlot

# Atom
Ip, d = hydrogen

# Field
λ = 800e-9 # m
I = 1e14 # W/cm²
ndt = 200
ω,T,A,E = field(λ,I)
tlims = (1,5,0)

I_au = I/3.5094452e16

function cutoff(Ip)
    Up = I_au/(4ω^2)
    c = (3.17Up + Ip)/ω
    a = gca()
    a2 = a[:twiny]()
    a2[:set_xlim](a[:get_xlim]())
    a2[:set_xticks]([c])
    a2[:set_xticklabels](["Cut-off"])
    sca(a)
end

@time x,t = propagate(A,E,Ip,d,tlims,T,ndt)
freq = fftshift(fftfreq(length(t),ndt))
X = fftshift(fft(hanning(length(t)).*x))

At = A(t)
Et = E(t)

figure("SFA HHG")
clf()
subplot(411)
plot(t/T, At, label=L"$A(t)$")
plot(t/T, Et, label=L"$E(t)$")
legend(framealpha=0.75)
subplot(412)
plot(t/T, x)
ylabel(L"$x(t)$")
xlabel(L"$t/T$")
subplot(413)
semilogy(freq,abs(X))
ylabel(L"$|x(\omega)|$")
xlim(0,50)
cutoff(Ip)
subplot(414)
plot(freq,unwrap(angle(X)))
xlabel(L"$\omega/\omega_0$")
ylabel(L"$\arg\{x(\omega)\}$")
xlim(0,50)
cutoff(Ip)
tight_layout()
