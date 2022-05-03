module ODE_mio

# Libraries this package relies upon
using FFTW, ProgressMeter

#The objects 'exported' here become available for the user after importing this package
export integrate_ode, derivative, rungekutta4, euler

# Its awkward to have everything in the same file.
# Lets split the package into multiple files!
include("include/integrator.jl")
include("include/derrivative.jl")

end #module

