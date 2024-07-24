#=
Samuel Nasman (BYU FLOW Lab)
Optimizing a wing to minimize induced drag
17 July 2024
=#

using VortexLattice, SNOW

include("./VLM_parameters.jl")

function min_C_D!(g, x)
    # objective
    grid, surface = wing_to_surface_panels(-0.25*x, yle, zle, x, theta, phi, ns, nc; fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=mirror)
    surfaces = [surface]
    Sref = ref_area(yle, x)
    cref = Sref / bref
    ref = Reference(Sref, cref, bref, rref, Vinf)
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
    CF, CM = body_forces(system; frame=Wind())
    CD, CY, CL = CF
    L = CL * 0.5 * Vinf^2 * Sref

    f = far_field_drag(system) * 0.5 * Vinf^2 * Sref

    # constraints
    g[1] = L
    g[2:(ns+1)] = diff(x)

    # chord_data[iter, :] = transpose(x)
    # lift_data[iter] = L
    # S_ref_data[iter] = Sref
    # drag_data[iter] = CD
    # drag_data_iff[iter] = f
    # global iter += 1

    return f
end

x0 = 10.0 * ones(ns+1) # starting point
lx = 0.00001 * ones(ns+1) # lower bounds on x
ux = Inf * ones(ns+1) # upper bounds on x
ng = ns+1  # number of constraints
lg = vcat(1.7, -Inf*ones(ns))  # lower bounds on g
ug = vcat(Inf, -0.000001*ones(ns)) # Inf * ones(ns+1)  # upper bounds on g
ip_options = Dict(
    # "max_iter" => 3000
)
options = Options(solver=IPOPT(ip_options))  # choosing IPOPT solver

# iter = 1
# max_iter = 20000
# chord_data = zeros(max_iter, ns+1)
# lift_data = zeros(max_iter)
# S_ref_data = zeros(max_iter)
# drag_data = zeros(max_iter)
# drag_data_iff = zeros(max_iter)


xopt, fopt, info = minimize(min_C_D!, x0, ng, lx, ux, lg, ug, options)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)

grid, surface = wing_to_surface_panels(-0.25*xopt, yle, zle, xopt, theta, phi, ns, nc; fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=mirror)
surfaces = [surface]
Sref = ref_area(-0.25*xopt, xopt)
cref = Sref / bref
rref = [0.0, 0.0, 0.0]
ref = Reference(Sref, cref, bref, rref, Vinf)
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
CF, CM = body_forces(system; frame=Wind())
CD, CY, CL = CF
CDiff = far_field_drag(system)
L = 0.5 * Vinf^2 * CL * Sref
write_vtk("test", surfaces; symmetric)