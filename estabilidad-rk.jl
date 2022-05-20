using SpecialFunctions, PyPlot # Importamos las funciones a utilizar 

# Definimos las funciones
function g2(z)
	1 .+ z .+ z.^2/2
end

function g4(z)
	1 .+ z .+ z.^2/2 .+ z^3/6 .+ z^4/24
end

# Definimos el eje x
x = -3:0.1:3;

# Definimos z y los valores absolutos de g
y = x;
z = x .+ im*y';
Zr2 = @. abs(g2(z));
Zr4 = @. abs(g4(z));

# Creamos el grafico
fig, ax = subplots()
ax.contour(x, y, Zr2', levels=[1], colors="blue")
ax.contour(x, y, Zr4', levels=[1], colors="red")
annotate("RK2", xy=[-2.4;0])
annotate("RK4", xy=[0.3;2])
axis([-3,3,-3,3])
title("Regiones de estabilidad")
xlabel("Re(z)")
ylabel("Im(z)")
savefig("fig-estabilidad-RK.pdf")
