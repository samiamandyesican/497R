using VortexLattice, Plots, Plots.PlotMeasures
localpath = @__DIR__
cd(localpath)
include("vortex_lattice_functions_wing_only.jl")

#= 
To analyze the effect of aspect ratio on inviscid span efficiency, I'm going to
sweep across a range of spans (by adjusting 'yle') while keeping everything else constant. 
=#


# Define Panels ------------------------------------------------------------
# Three currently available discretization scheme options: Uniform(), Sine() and Cosine(). 
ns = 20 # number of spanwise panels
nc = 8  # number of chordwise panels
spacing_s = Sine() # spanwise discretization scheme
spacing_c = Uniform() # chordwise discretization scheme

# Define freestream -------------------------------------------------------
alpha = 1.0*pi/180 # angle of attack
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
# we can model other half of wing using 'symmetric = true', but it would give incorect lateral stability, so we instead already mirrored our geometry across X-Z plane
symmetric = false

# Define Geometry --------------------------------------------------------
xle = [0.0, 0.4] # leading edge x-position
zle = [0.0, 0.0] # leading edge z-position
chord = [2.2, 1.8] # chord length
theta = [2.0*pi/180, 2.0*pi/180] # twist (in radians)
phi = [0.0, 0.0] # section rotation about the x-axis
fc = fill((xc) -> 0, 2) # camberline function for each section (y/c = f(x/c))

# Define constant reference geometry --------------------------------------
cref = ref_chord(chord)  # reference chord
rref = center_of_g(chord) # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity 
fs = Freestream(Vinf, alpha, beta, Omega)


# Define iteration parameters and allocate matrices -------------------------------------------------
min_half_span = 0.5
max_half_span = 100
reiterations = 20
step = (max_half_span - min_half_span) / (reiterations)
AR = zeros(reiterations+1) # aspect ratios 
eff = zeros(reiterations+1) # inviscid span efficiencies
C_Di = zeros(reiterations+1) # induced drag coefficient
C_D = zeros(reiterations+1) # drag coefficient
C_L = zeros(reiterations+1) # lift coefficient
lift_drag_ratio = zeros(reiterations+1) # lift to drag coefficient ratio
L = zeros(reiterations+1) # lift
D = zeros(reiterations+1) # drag


for i in 0:reiterations

    yle = [0.0, min_half_span + i * step] # leading edge y-position (Note that we will sweep through a range of these)
    # Generate Grid / Lifting Surface ------------------------------------------
    # Note that the grid is not needed for analysis, just there for convenience.
    # Also note that we could generate a lifting surface from a pre-existing grid using 'grid_to_surface_panels'.
    # 'mirror' mirrors geometry across the X-Y plane.
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc; fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)
    # Combine all surfaces in a single vector (in this case we only have one surface)
    surfaces = [surface]

    # Define reference parameters ----------------------------------------------
    Sref = ref_area_trap(yle, chord)
    bref = ref_span(yle) # reference span
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # Steady state analysis ----------------------------------------------------
    # performs near-field analysis on each panel unless kwarg 'near_field_analysis' is changed
    # Determines derivatives with respect to freestream variables unless specified through kwarg 'derivatives'
    # kwarg 'symmetric' is not strictly necessary, since by default it is set to false for each surface.
    system = steady_analysis(surfaces, ref, fs; symmetric)

    # Extract force/moment coefficients -----------------------------------------
    # returned in the reference fram specified by kwarg 'frame', with the body being the default reference
    # only returns sensible results if near field analysis was performed on 'system'
    CF, CM = body_forces(system; frame=Wind())
    # extract aerodynamic forces
    CD, CY, CL = CF
    Cl, Cm, Cn = CM
    # Numerical noise often corrupts drag estimates so it is often more accurate to compute drag in the farfield on the Trefftz plane
    CDiff = far_field_drag(system)

    AR[i+1] = aspect_ratio(yle, chord)
    eff[i+1] = wing_efficiency(CL, CDiff, yle, chord)
    C_Di[i+1] = CDiff
    C_D[i+1] = CD
    C_L[i+1] = CL
    lift_drag_ratio[i+1] = CL / CD
    L[i+1] = force(CL, Sref, Vinf)
    D[i+1] = force(CD, Sref, Vinf)

    if i == 0
        write_vtk("low-AR", surfaces; symmetric)
    elseif i == reiterations
        write_vtk("high-AR", surfaces; symmetric)
    end

end

plt1 = plot(AR, eff; xlabel="Aspect Ratio", ylabel="Inviscid Span Efficiency", leg=false)
plt2 = plot(AR, C_Di; xlabel="Aspect Ratio", ylabel="Induced Drag Coefficient", leg=false)
plt3 = plot(AR, C_D; xlabel="Aspect Ratio", ylabel="Drag Coefficient", leg=false)
plt4 = plot(AR, C_L; xlabel="Aspect Ratio", ylabel="Lift Coefficient", leg=false)
plt5 = plot(AR, lift_drag_ratio; xlabel="Aspect Ratio", ylabel="Lift / Drag Coefficient Ratio", leg=false)
plt6 = plot(AR, L; xlabel="Aspect Ratio", ylabel="Lift", leg=false)
plt7 = plot(AR, D; xlabel="Aspect Ratio", ylabel="Drag", leg=false)

plt = plot(plt1, plt2, plt3, plt4, plt5, plt6, plt7; layout=(7,1))
display(plt)

# plt = plot(AR, [eff C_Di C_L lift_drag_ratio L D], xlabel="Aspect Ratio", markershape=:x, label=["Inviscid Span Efficiency" "Induced Drag Coefficient" "Lift Coefficient" "Lift / Drag Coefficient Ratio" "Lift" "Drag"], layout=(6,1), leg=:best, size=(700,2500), left_margin=20.0mm)
# display(plt)
# savefig("aspect_ratio-test.png")

nothing