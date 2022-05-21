using Plots, ProgressMeter

function derivative(fx, x; order=1)
	Nx = length(x)
	dx = x[2]-x[1]

	ik = im .* rfftfreq(Nx, 2pi/dx)

	return irfft(rfft(fx) .* ik.^order, Nx)
end

function rungekutta4(dydt, y, t, dt)
	k1 = dt .* dydt(y, t)
	k2 = dt .* dydt(y .+ k1./2, t + dt/2)
	k3 = dt .* dydt(y .+ k2./2, t + dt/2)
	k4 = dt .* dydt(y .+ k3   , t + dt)

	return @. y + (k1 + 2k2 + 2k3 + k4) / 6 # y(t + dt)
end

function integrate_ode(dydt, y0, Nt, dt)

	# Time axis
	time = range(dt, step=dt, length=Nt)

	# Define an accumulator for the i-th value at the j-th time
	trajectory = Array{Float64}(undef, length(y0), Nt)

	y = copy(y0)
	progressbar = Progress(Nt, desc="Integrating...")
	for (i, t) in enumerate(time)

			# Obtain y(t+dt) from y(t) and dydt(y, t)
			y = rungekutta4(dydt, y, t, dt)

			# Save the i-th value of y, y(t=i*dt)
			trajectory[:, i] .= y

			# Update progressbar
			next!(progressbar)
	end
	return time, trajectory
end

# Define the space
Nx = 256
xmin, xmax = 0.0, 2pi
x = range(xmin, step=(xmax-xmin)/Nx, length=Nx)

# Temporal discretization
dt = 3e-3
Nt = round(Int, 1.5/dt)

# Initial conditions
u = sin.(x)

# Differential equation
dudt(u, t) = -0.5derivative(u.^2, x) + 1e-3*derivative(u, x, order=2)

# Integrate initial conditions with ODE.integrate_ode
time, ut = integrate_ode(dudt, u, Nt, dt)

# Make animation with only some instants
instants_to_plot = 1:Nt

progressbar = Progress(length(instants_to_plot), desc="Animating...")
anim = @animate for i in instants_to_plot
    p = plot()

    # Plot initial conditions
    plot!(p, x, ut[:, 1], l = (:gray, :dash), label="\$u(x, 0)\$")
    # Plot ut at current time
    plot!(p, x, ut[:, i], l = (:black), label="\$u(x, t)\$")
    # Add other labels
    plot!(p,
          xlabel = "\$ x\$",
          ylabel = "\$ u\$",
          title = "t = $(round(time[i], digits=2))",
          )
end

# Save animation as mp4
mp4(anim, "burgers.mp4")
