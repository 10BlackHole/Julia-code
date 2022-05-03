using ODE_mio
using Plots, ProgressMeter

#Definde the space
Nx = 256
xmin, xmax = 0.0, 2pi
x = range(xmin, step=(xmax-xmin)/Nx, length=Nx)

#Temporal discretization
Nt = 2000
dt = 1e-3

#Initial conditions
u = sin.(x)

# Differential equation
v(x) = 1.0
dudt(u, t)= -v(x) .* derivative(u, x)

# Integrate initial conditions with ODE,integrate_ode
time, ut = integrate_ode(dudt, u, Nt, dt)

# Make animation with only some instants
instants_to_plot = 1:10:Nt

progressbar = Progress(length(instants_to_plot), desc="Animating...")
anim = @animate for i in instants_to_plot
    p = plot()

    # Plot initial conditions
    plot!(p, x, ut[:, 1], l=(:gray, :dash), label="\$u(x,0)\$")
    # Plot it at current time
    plot!(p, x, ut[:, i], l = (:black), label="\$u(x,t)\$")
    # Add other labels
    plot!(p,
         xlabel =  "\$x\$",
         ylabel = "\$u\$",
         title = "t = $(round(time[i],digits=2))",
        )
end
