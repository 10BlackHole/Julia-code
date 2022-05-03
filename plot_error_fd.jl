using PyPlot

function f_deriv(f,h,x) #forward difference
	deriv = (f(x+h)-f(x))/(h)
	return deriv
end

function b_deriv(f,h,x) # backward difference
	deriv = (f(x)-f(x-h))/h
	return deriv
end

function c_deriv(f,h,x) # central difference
	deriv = (f(x+h)-f(x-h))/(2h)
	return deriv
end

f(x) = cos(2pi*x) + sin(4pi*x)
df(x) = -2pi*sin(2pi*x) + 4pi*cos(4pi*x)

h = 1e-2
x0 = 0.5

err_f = abs(f_deriv(f,h,x0) - df(x0))
err_b = abs(b_deriv(f,h,x0) - df(x0))
err_c = abs(c_deriv(f,h,x0) - df(x0))

# println(err_f)
# println(err_b)
# println(err_c)

np = 50
hx = 10 .^ range(-10, -2, length = np)

error_cp = zeros(np)
error_fp = similar(error_cp)
error_bp = similar(error_cp)

for i=1:np
	h = hx[i]
	error_fp[i] = abs(f_deriv(f,h,x0) - df(x0))
	error_bp[i] = abs(b_deriv(f,h,x0) - df(x0))
	error_cp[i] = abs(c_deriv(f,h,x0) - df(x0))
end

fig, ax = subplots()
ax.loglog(hx, error_fp, "o", label=L"adelantada", alpha=0.5)
ax.loglog(hx, error_bp, "-", label=L"retrasada")
ax.loglog(hx, error_cp, ".", label=L"centrada") 
ax.legend(loc="upper left")
axis([1e-10, 1e-2, 1e-15, 1e0])
title("Error de la 1era derivada con diferencias finitas")
ylabel(L"|f'_{FD} - f'_{exact}|")
xlabel(L"h")
# savefig("Error_FD.pdf")
