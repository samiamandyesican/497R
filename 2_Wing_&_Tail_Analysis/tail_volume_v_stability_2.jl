using VortexLattice
include("vortex_lattice_functions_tail.jl")

#=
I'm going to sweep across a range of tailplane geometries to examine the effect of tail volume ratio on stability
=#

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 0.0]
chord = [1.0, 0.6]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 12
nc = 6
spacing_s = Sine()
spacing_c = Uniform()
mirror = true

# reference geometry
Sref = ref_area_trap(yle, chord)
cref = ref_chord(chord)
bref = ref_span(yle)
rref = [cref, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream
alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

symmetric = [false, false, false]

trans = [10.0, 0.0, 0.0]

# define iteration parameters
iterations = 20 # number of iterations
min_size_v = 0.5 # magnification of vertical stabilizer
max_size_v = 2.5 # magnification of vertical stabilizer
min_size_h = 1.0 # magnification of horizontal stabilizer
max_size_h = 1.0 # magnification of horizontal stabilizer

    # horizontal stabilizer initial values
    xle_h = [0.0, 0.14]
    yle_h = [0.0, 2.0]
    zle_h = [0.0, 0.0]
    chord_h = [0.7, 0.42]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    fc_h = fill((xc) -> 0, 2) #camberline function for each section
    ns_h = 6
    nc_h = 3
    spacing_s_h = Sine()
    spacing_c_h = Uniform()
    mirror_h = true

    xle_h_steps = get_steps_vec(min_size_h * xle_h, max_size_h * xle_h, iterations)
    yle_h_steps = get_steps_vec(min_size_h * yle_h, max_size_h * yle_h, iterations)
    chord_h_steps = get_steps_vec(min_size_h * chord_h, max_size_h * chord_h, iterations)
    ns_h_steps = get_steps_num(min_size_h * ns_h, max_size_h * ns_h, iterations)
    nc_h_steps = get_steps_num(min_size_h * nc_h, max_size_h * nc_h, iterations)

    # vertical stabilizer initial values
    xle_v = [0.0, 0.14]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, 2.0]
    chord_v = [0.7, 0.42]
    theta_v = [0.0, 0.0]
    phi_v = [0.0, 0.0]
    fc_v = fill((xc) -> 0, 2) #camberline function for each section
    ns_v = 5
    nc_v = 3
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

    xle_v_steps = get_steps_vec(min_size_v * xle_v, max_size_v * xle_v, iterations)
    zle_v_steps = get_steps_vec(min_size_v * zle_v, max_size_v * zle_v, iterations)
    chord_v_steps = get_steps_vec(min_size_v * chord_v, max_size_v * chord_v, iterations)
    ns_v_steps = get_steps_num(min_size_v * ns_v, max_size_v * ns_v, iterations)
    nc_v_steps = get_steps_num(min_size_v * nc_v, max_size_v * nc_v, iterations)

# initialize matrices for data

labels = ["V_th" "V_tv" "CDa" "CYa" "CLa" "Cla_" "Cma" "Cna" "CDb" "CYb" "CLb" "Clb_" "Cmb" "Cnb" "CDp" "CYp" "CLp" "Clp_" "Cmp" "Cnp" "CDq" "CYq" "CLq" "Clq_" "Cmq" "Cnq" "CDr" "CYr" "CLr" "Clr_" "Cmr" "Cnr"]
data = zeros(iterations,32)

# generate surface panels for wing
wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc; mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)


for i in 1:iterations

    # horizontal stabilizer
    xle_h = xle_h_steps[:, i]
    yle_h = yle_h_steps[:, i]
    chord_h = chord_h_steps[:, i]
    ns_h = ns_h_steps[i]
    nc_h = nc_h_steps[i]

    # vertical stabilizer
    xle_h = xle_v_steps[:, i]
    zle_v = zle_v_steps[:, i]
    chord_v = chord_v_steps[:, i]
    ns_v = ns_v_steps[i]
    nc_v = nc_v_steps[i]

    # generate surface panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(hgrid, trans)
    translate!(htail, trans)

    # generate surface panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vgrid, trans)
    translate!(vtail, trans)

    grids = [wgrid, hgrid, vgrid]
    surfaces = [wing, htail, vtail]
    surface_id = [1, 2, 3]

    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

    CF, CM = body_forces(system; frame=Wind())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    properties = get_surface_properties(system)

    # Extract stability derivatives -------------------------------------------------
    # Note that for sensible results, near-field analysis must have been performed
    dCFs, dCMs = stability_derivatives(system)
    # traditional names for each stability derivative
    CDa, CYa, CLa = dCFs.alpha
    Cla, Cma, Cna = dCMs.alpha
    CDb, CYb, CLb = dCFs.beta
    Clb, Cmb, Cnb = dCMs.beta
    CDp, CYp, CLp = dCFs.p
    Clp, Cmp, Cnp = dCMs.p
    CDq, CYq, CLq = dCFs.q
    Clq, Cmq, Cnq = dCMs.q
    CDr, CYr, CLr = dCFs.r
    Clr, Cmr, Cnr = dCMs.r
    V_th = tail_volume_h(yle, yle_h, chord, chord_h, trans, rref)
    V_tv = tail_volume_v(yle, xle_v, zle_v, chord, chord_v, trans, rref)

    data[i, :] = [V_th V_tv CDa CYa CLa Cla Cma Cna CDb CYb CLb Clb Cmb Cnb CDp CYp CLp Clp Cmp Cnp CDq CYq CLq Clq Cmq Cnq CDr CYr CLr Clr Cmr Cnr]

    # if i ==1
    #     write_vtk("wing-tail-small-h", surfaces, properties; symmetric)
    # elseif i == iterations
    #     write_vtk("wing-tail-large-h", surfaces, properties; symmetric)
    #     global system = system
    # end

end

data_l = vcat(labels, data)

include("stability_plots.jl")
generate_plots(labels, data, "vertical_scaling_only")

nothing