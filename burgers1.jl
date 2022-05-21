using FFTW
using PyPlot

function wavevector(L, Nx);
	knx = collect(1:Nx) .- (Nx/2 + 1) # integer wavenumbers for X-space
	kx = 2pi*fftshift(knx)/L # space-wavenumber vector
end

function rg4(f, y0, h, kx, D)
	n = length(y0)
	y = zeros(n)
		k1 = f(y0, kx, D)
		k2 = f(y0 .+ k1 * h/2, kx, D)
		k3 = f(y0 .+ k2 * h/2, kx, D)
		k4 = f(y0 .+ k3 * h, kx, D)
		y = y0 .+ (h/6) * (k1 + 2*k2 + 2*k3 + k4)
	return y
end

function F(y, kx, D)
	du2x, d2ux = deriv(y, kx)
	lder = -0.5*du2x .+ D*d2ux
end

function deriv(u, kx)
	uk = fft(u)
	u2k = fft(u.^2)
	du2x = real(ifft(im*kx.*u2k))
	d2ux = real(ifft(-kx.^2 .*uk))
	return du2x, d2ux
end

Nx = 128 # 
L = 2pi # 
dt = 0.003
h = L/Nx #as celdas
x = (1:Nx)*h # la variable independiente espacial x
u0 = sin.(x) 
kx = wavevector(L, Nx) 
D = 0.01
Tfin = 1.5
u1 = u0 
u = similar(u0)
for i in 0:dt:Tfin
	u = rg4(F, u1, dt, kx, D)
	u1 = u
end

plot(x, u0)
plot(x, u)
