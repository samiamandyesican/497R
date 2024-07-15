using VortexLattice, Plots
localpath = @__DIR__
cd(localpath)
include("vortex_lattice_functions_wing_only.jl")

# Define Panels ------------------------------------------------------------
ns = 20 # number of spanwise panels
nc = 8  # number of chordwise panels

# Define freestream -------------------------------------------------------
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
# we can model other half of wing using 'symmetric = true', but it would give incorect lateral stability, so we instead already mirrored our geometry across X-Z plane
symmetric = false

# Define Geometry --------------------------------------------------------
xle = [0.0, 0.4] # leading edge x-position
yle = [0.0, 7.5] # leading edge y-position
zle = [0.0, 0.0] # leading edge z-position
chord = [2.2, 1.8] # chord length [root, tip]
theta = [2.0*pi/180, 2.0*pi/180] # twist (in radians)
phi = [0.0, 0.0] # section rotation about the x-axis
fc = fill((xc) -> 0, 2) # camberline function for each section (y/c = f(x/c))

# Define constant reference geometry --------------------------------------
cref = ref_chord(chord)  # reference chord
rref = center_of_g(chord) # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity 

# leading edge y-position (Note that we will sweep through a range of these)
# Generate Grid / Lifting Surface ------------------------------------------
# Note that the grid is not needed for analysis, just there for convenience.
# Also note that we could generate a lifting surface from a pre-existing grid using 'grid_to_surface_panels'.
# 'mirror' mirrors geometry across the X-Y plane.
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc; fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)
# Combine all surfaces in a single vector (in this case we only have one surface)
surfaces = [surface]

# Define reference parameters ----------------------------------------------
Sref = ref_area_trap(yle, chord) # reference area
bref = ref_span(yle) # average (reference) span
ref = Reference(Sref, cref, bref, rref, Vinf)


# Define iteration parameters and allocate matrices -------------------------------------------------
iterations = 20 # number of iterations
alpha_min = 0.0*pi/180 # min angle of attack
alpha_max = 360.0*pi/180 # max angle of attack
alpha_step = get_step_float(alpha_min, alpha_max, iterations) # create a vector of 20 equally spaced alpha angles from max to min.
data = zeros(iterations, 2) # allocate matrix to hold data

for i in 1:iterations

    fs = Freestream(Vinf, alpha_step[i], beta, Omega) # define freestream using alpha

    # Steady state analysis ----------------------------------------------------
    # performs near-field analysis on each panel unless kwarg 'near_field_analysis' is changed
    # Determines derivatives with respect to freestream variables unless specified through kwarg 'derivatives'
    # kwarg 'symmetric' is not strictly necessary, since by default it is set to false for each surface.
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

    # Extract force/moment coefficients -----------------------------------------
    # returned in the reference frame specified by kwarg 'frame', with the body being the default reference
    # only returns sensible results if near field analysis was performed on 'system'
    CF, CM = body_forces(system; frame=Wind())

    # extract aerodynamic forces
    CD, CY, CL = CF
    Cl, Cm, Cn = CM
    # Numerical noise often corrupts drag estimates so it is often more accurate to compute drag in the farfield on the Trefftz plane
    CDiff = far_field_drag(system)

    # save data to data matrix
    data[i, 1] = alpha_step[i] * 180 / pi
    data[i, 2] = CL

end

# plot data and save results (replacing file of same name if it exists)
plt = plot(data[:,1], data[:,2], xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", leg=:best)
save_path = "alpha_v_CL"
rm(save_path, force=true)
savefig(plt, save_path)

nothing