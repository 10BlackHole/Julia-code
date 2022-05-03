# Integrators should take the ODE 'dydt'
# an array 'y' representign the function at time 't'
# the time 't',  and the time step 'dt'
# They will return an array of the function at time (t + dt)

function rungekutta4(dydt, y, t, dt)
    k1 = dt .* dydt(y, t)
    k2 = dt .* dydt(y .+ k1./2, t + dt/2)
    k3 = dt .* dydt(y .+ k2./2, t + dt/2)
    k4 = dt. * dydt(y .+ k3, t + dt)
    return @. y + (k1 + 2k2 + 2k3 + k4) / 6 #y(t + dt)
end

#----------------------------------------------------------
#First order Euler integrator

function euler(dydt, y, t, dt)
    return @. y + dt * dydt(y, t) #y(t + dt)
end


#-------------------------------------------------------------
# Integrates de ODE 'dydt', startign from the array with initial conditions 'y0' during 'Nt'
# iterations with a time step 'dt', using the integrator 'integrator'.
# Returns a tuple '(time, trajectory)', where 'time' is the time axis ranging from 'dt' to 'dt*Nt' with steps 'dt',
# and  'trajectory' is an accumulator, 'trajectory[i, j]' es the 'i-tj' of 'y' at the 'j-th' time.
function integrate_ode(dydt, y0, Nt, dt; integrator=euler)

    # Time axis
    time = range(dt, step=dt, length=Nt)

    # Define an accumulator for the i-th value at the j-th time
    trajectory = Array{Float64}(undef, length(y0), nT)

    y = copy(y0)
    progressbar = Progress(Nt, desc="Integrating...")
    for (i, t) in enumerate(time)

        # Obtain y(t +dt) from y(t) and dydt(y, t)
        y = integrator(dydt, y, t, dt)

        #Save the i-th value of y, y(t = i*dt)
        trajectory[:, i]. = y

        # Update progressbar
        next!(progressbar)
    end

    return time, trajectory
end

