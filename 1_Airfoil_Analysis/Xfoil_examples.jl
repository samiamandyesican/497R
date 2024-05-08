#=
FLOW Lab 497R
Samuel Nasman
4 May 2024

# Examples from https://flow.byu.edu/Xfoil.jl/stable/guide/#guide
# =#

using Xfoil, Plots, Printf

pyplot()


# ------------------ LOADING AIRFOIL GEOMETRY ----------------------- #
# read airfoil coordinates from a file
function read_coordinates(filename)
    open(filename, "r") do f
        x = Float64[]
        y = Float64[]
        i = 0
        
        for line in eachline(f)
            # skip over title
            if i in 0:0
                i += 1
            
            # read entries
            else
                entries = split(chomp(line))
                push!(x, parse(Float64, entries[1]))
                push!(y, parse(Float64, entries[2]))
            end
        end

        return x, y    

    end
end


function load_coordinates(x, y)
    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x, y)
end


function plot_airfoil_geometry(x, y)
    # plot the airfoil geometry
    scatter(x, y, label="", framestyle=:none, aspect_ratio=1.0, show=true)
end


"""
    insert_airfoil_geometry(filename)

Load coordinates of an airfoil from a file into XFoil and then plot the coordinates. 
"""
function insert_airfoil_geometry(filename)
    x, y = read_coordinates(filename)
    load_coordinates(x, y)
    plot_airfoil_geometry(x, y)
end

x, y = read_coordinates("naca2412.txt")
insert_airfoil_geometry("naca2412.txt")


#----------------REFINING AIRFOIL GEOMETRY---------------------#
"""
    pane(; kwargs...)

Run XFOIL's PANE command (Repanel the input airfoil)

# Arguments
- `npan::Integer=140`: Number of panel nodes
- `cvpar::Float64=1.0`: Panel bunching parameter
- `cterat::Float64=0.15`: TE/LE panel density ratio
- `ctrrat::Float64=0.2`: Refined-area/LE panel density ratio
- `xsref1::Float64=1.0`: Top side refined area x/c limits
- `xsref2::Float64=1.0`:
- `xpref1::Float64=1.0`: Bottom side refined area x/c limits
- `xpref2::Float64=1.0`:
"""
# repanel using XFOIL's `PANE` command
xr, yr = Xfoil.pane()

# plot the refined airfoil geometry
scatter(xr, yr, label="", framestyle=:none, aspect_ratio=1.0, overwrite_figure=false, show=true)


#----------------------DEFINING OPERATING CONDITIONS----------------------#
"""
For inviscid analyses, only the angle of attack is required. For viscous analyses, the Reynolds number must also be specified.
"""
# set operating conditions
alpha = -9:1:14 # range of angle of attacks, in degrees
re = 1e5 # Reynolds number


#--------------------------AIRFOIL ANALYSIS----------------------------------#
"""
The solve_alpha function may now be used to perform an analysis to obtain the airfoil coefficients c_l, c_d, c_d_p, and c_m.
Note that c_d_p is profile drag.
Skin friction drag may be obtained by c_d_f = c_d - c_d_p.

Note that the order in which viscous analyses are performed matters since XFOIL uses boundary layer parameters 
    corresponding to the last previously converged solution as an initial guess 
    when solving for the current boundary layer parameters. 
This behavior can be disabled by passing the keyword argument pair reinit=true to solve_alpha
"""
# initialize outputs
n_a = length(alpha)
c_l = zeros(n_a)
c_d = zeros(n_a)
c_dp = zeros(n_a)
c_m = zeros(n_a)
converged = zeros(Bool, n_a)

# determine airfoil coefficients across a range of angle of attacks
for i = 1:n_a
    c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
end

# print results
println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
for i = 1:n_a
    @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
end

# plot results
plot(alpha, c_l, label="", xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", 
    overwrite_figure=false, show=true)
plot(alpha, c_d, label="", xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient",
    overwrite_figure=false, show=true)
plot(alpha, c_m, label="", xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient",
    overwrite_figure=false, show=true)


#-----------------------SENSITIVITY ANALYSIS---------------------------------#
"""
Suppose we want to find the derivative of c_l, c_d, and c_m with respect to the angle of attack. 
One approach to calculate these derivatives would be to use the finite difference method.
"""

# set operating conditions
alpha = -9:1:14
re = 1e5
mach = 0.0

# set step size
h = 1e-6

# initialize outputs
n_a = length(alpha)
c_l_a = zeros(n_a)
c_d_a = zeros(n_a)
c_dp_a = zeros(n_a)
c_m_a = zeros(n_a)
converged = zeros(Bool, n_a)

for i = 1:n_a
    c_l1, c_d1, c_dp1, c_m1, converged[i] = Xfoil.solve_alpha(alpha[i], re; mach, iter=100, reinit=true)
    c_l2, c_d2, c_dp2, c_m2, converged[i] = Xfoil.solve_alpha(alpha[i]+h, re; mach, iter=100, reinit=true)
    c_l_a[i] = (c_l2 - c_l1)/h * 180/pi
    c_d_a[i] = (c_d2 - c_d1)/h * 180/pi
    c_m_a[i] = (c_m2 - c_m1)/h * 180/pi
end

# print results
println("Angle\t\tdClda\t\tdCdda\t\tdCmda\t\tConverged")
for i = 1:n_a
  @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l_a[i],c_d_a[i],c_m_a[i],converged[i])
end

"""
A better approach might be to use the complex step method. 
    To use this approach, however, we must use the complex-step enabled version of XFOIL provided by this package.
The complex-step version of each function is denoted by appending _cs to each function name. 
    Note that there is no interaction between the two versions of XFOIL wrapped by this package, so if you wish to use the complex step version of the code you must append _cs to all function names.
For the complex step method to work, we also need to ensure that the imaginary portion of variables from previous 
    iterations does not affect the solution for the current iteration. This may be achieved by setting reinit=true when calling solve_alpha.

        Note that since XFOIL was not originally designed for sensitivity analysis, there is a distinct possibility that sensitivities may be non-physical. 
            We therefore always recommend checking that computed sensitivities are realistic. 
            Sometimes adjusting the number of panels or the step size will fix the computed sensitivies, 
            but if all else fails a coarse step size (e.g. 1.0 degree) may be used with central finite differencing to artificially smooth 
                the computed sensitivities.
"""

# read airfoil into XFOIL
Xfoil.set_coordinates_cs(x,y)

# repanel using XFOIL's `PANE` command
Xfoil.pane_cs()

# set operating conditions
alpha = -9:1:14
re = 1e5
mach = 0.0

# set step size
h = 1e-12im

# initialize outputs
n_a = length(alpha)
c_l_a = zeros(n_a)
c_d_a = zeros(n_a)
c_dp_a = zeros(n_a)
c_m_a = zeros(n_a)
converged = zeros(Bool, n_a)

for i = 1:n_a
    c_l, c_d, c_dp, c_m, converged[i] = Xfoil.solve_alpha_cs(alpha[i]+h, re; mach, iter=100, reinit=true)
    c_l_a[i] = imag(c_l)/imag(h) * 180/pi
    c_d_a[i] = imag(c_d)/imag(h) * 180/pi
    c_m_a[i] = imag(c_m)/imag(h) * 180/pi
end

# print results
println("Angle\t\tdClda\t\tdCdda\t\tdCmda\t\tConverged")
for i = 1:n_a
  @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l_a[i],c_d_a[i],c_m_a[i],converged[i])
end



#-------------------------AUTOMATED ANGLE OF ATTACK SWEEP--------------------------#
"""
For performing angle of attack sweeps, the function alpha_sweep may also be used.
"""

# set operating conditions
alpha = -9:1:14
re = 1e5

c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(x, y, alpha, re, iter=100, zeroinit=false, printdata=true, reinit=true)

"""
A version of alpha_sweep has also been implemented for use with the complex step version of XFOIL.
"""

# set operating conditions
alpha = -9:1:14
re = 1e5
mach = 0.0

# set step size
h = 1e-20im

c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep_cs(x, y, alpha .+ h,
    re, mach=mach, iter=100, zeroinit=false, printdata=true, reinit=true)

println("Angle\t\tdClda\t\tdCdda\t\tdCmda\t\tConverged")
for i = 1:length(alpha)
    @printf("%8f\t%8f\t%8f\t%8f\t%d\n", alpha[i], imag(c_l[i])/imag(h)*180/pi, imag(c_d[i])/imag(h)*180/pi, imag(c_m[i])/imag(h)*180/pi, converged[i])
end