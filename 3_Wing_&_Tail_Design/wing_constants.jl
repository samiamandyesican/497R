include("./ref_geometry_functions.jl")

# Define Panels ------------------------------------------------------------
spacing_s = Uniform()
spacing_c = Uniform()
#= we can model other half of wing using 'symmetric = true', but it would give 
incorect lateral stability, so we instead will mirror our geometry across X-Z plane 
using 'mirror=true' =#
symmetric = false

# Define Chord Independent Geometry  --------------------------------------------------------
function define_geometry(ns)
    yle = stepped_vec(0.0, 4.0, ns+1) # leading edge y-positions
    zle = zeros(Float64, ns+1) # leading edge z-positions
    theta = 0.0 * pi/180 * ones(ns+1) # twist (in radians)
    phi = zeros(Float64, ns+1) # section rotations about the x-axis
    fc = fill((xc) -> 0, ns+1) # camberline function for each section (y/c = f(x/c))

    # chord independent reference geometry
    bref = 2.0 * yle[ns+1] # span
    return yle, zle, theta, phi, fc, bref
end
rref = [0.0, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
mirror = true # 'mirror' mirrors geometry across the X-Y plane.

# Define freestream -------------------------------------------------------
alpha = 5.0 * pi/180 # angle of attack
beta = 0.0 * pi/180 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
Vinf = 1.0 # reference velocity 
fs = Freestream(Vinf, alpha, beta, Omega)