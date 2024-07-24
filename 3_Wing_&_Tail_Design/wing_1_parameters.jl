include("./ref_geometry_functions.jl")

# Define Panels ------------------------------------------------------------
ns = 6 # number of spanwise panels
nc = 8  # number of chordwise panels
spacing_s = Uniform()
spacing_c = Uniform()
#= we can model other half of wing using 'symmetric = true', but it would give 
incorect lateral stability, so we instead will mirror our geometry across X-Z plane 
using 'mirror=true' =#
symmetric = false

# Define Geometry  --------------------------------------------------------
# chord = 10.0 * ones(ns+1) # chord lengths
# xle = -0.25*chord # leading edge x-positions (to have constantly aligned quarter chords)
yle = stepped_vec(0.0, 4.0, ns+1) # leading edge y-positions
zle = zeros(Float64, ns+1) # leading edge z-positions
theta = 0.0 * pi/180 * ones(ns+1) # twist (in radians)
phi = zeros(Float64, ns+1) # section rotations about the x-axis
fc = fill((xc) -> 0, ns+1) # camberline function for each section (y/c = f(x/c))
mirror = true # 'mirror' mirrors geometry across the X-Y plane.

# Define reference geometry --------------------------------------
bref = 2.0 * yle[ns+1] # span
# Sref = ref_area(yle, chord) # reference area
# cref = Sref / bref  # reference chord (mean geometric chord)
rref = [0.0, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)

# Define freestream -------------------------------------------------------
alpha = 5.0 * pi/180
beta = 0.0 * pi/180 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
Vinf = 1.0 # reference velocity 

# structure data
# geometry = Geometry(xle, yle, zle, chord, theta, phi, fc, mirror)
# ref = Reference(Sref, cref, bref, rref, Vinf)
fs = Freestream(Vinf, alpha, beta, Omega)
# panels = Panel(ns, nc, spacing_s, spacing_c, symmetric)
# wing = WingDef(geometry, ref, fs, panels)