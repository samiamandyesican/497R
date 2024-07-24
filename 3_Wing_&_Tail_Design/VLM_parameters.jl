# generate data structures 
# struct Geometry
#     xle::Vector
#     yle::Vector
#     zle::Vector
#     chord::Vector
#     theta::Vector
#     phi::Vector
#     fc::Vector
#     mirror::Bool
# end

# struct Panel
#     ns::Int64
#     nc::Int64
#     spacing_s
#     spacing_c
#     symmetric::Bool
# end

# struct WingDef
#     geometry::Geometry
#     ref::Reference
#     fs::Freestream
#     panels::Panel
# end

include("wing_1_parameters.jl")

# (; geometry, ref, fs, panels) = wing
# (; xle, yle, zle, chord, theta, phi, fc, mirror) = geometry
# (; ns, nc, spacing_s, spacing_c, symmetric) = panels

# #Generate Grid / Lifting Surface ------------------------------------------
# grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc; fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=mirror)
# # Combine all surfaces in a single vector (in this case we only have one surface)
# surfaces = [surface]

# # Steady state analysis ----------------------------------------------------
# # performs near-field analysis on each panel unless kwarg 'near_field_analysis' is changed
# # Determines derivatives with respect to freestream variables unless specified through kwarg 'derivatives'
# # kwarg 'symmetric' is not strictly necessary, since by default it is set to false for each surface.
# system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# # Extract force/moment coefficients -----------------------------------------
# # returned in the reference frame specified by kwarg 'frame', with the body being the default reference
# # only returns sensible results if near field analysis was performed on 'system'
# CF, CM = body_forces(system; frame=Wind())

# # extract aerodynamic forces
# CD, CY, CL = CF
# Cl, Cm, Cn = CM
# # Numerical noise often corrupts drag estimates so it is often more accurate to compute drag in the farfield on the Trefftz plane
# CDiff = far_field_drag(system)

# write_vtk("test", surfaces; symmetric)

nothing