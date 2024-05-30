using DataFrames, CSV, Plots, Xfoil, Printf

localpath = @__DIR__
cd(localpath)

# ------------------ LOADING AIRFOIL GEOMETRY ----------------------- #
# read airfoil coordinates from a file
function read_coordinates_X(filename)
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


# ------------------ Plotting Experimental Data ----------------------- #
# read airfoil coordinates from a file
function read_coordinates(filename)
    open(filename, "r") do f
        global x = Float64[]
        global y = Float64[]
        
        for line in eachline(f)
            entries = split(chomp(line), ", ")
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))        
        end

        return x, y    
    
    end
end

#----------------------DEFINING OPERATING CONDITIONS----------------------#
"""
For inviscid analyses, only the angle of attack is required. For viscous analyses, the Reynolds number must also be specified.
"""
# set operating conditions
alpha = -20:1:25 # range of angle of attacks, in degrees
re = 5.7e6 # Reynolds number
mach = 0.13


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

# # print results
# println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
# for i = 1:n_a
#     @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
# end

# for re in [1e6, 1e7, 1e8, 1e9]
#     # determine airfoil coefficients across a range of angle of attacks
#     c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(xr, yr, alpha, re, iter=100, zeroinit=false, printdata=true, reinit=true)
#     plot!(plt, alpha, c_m, label="Re = $re")
# end

# # determine airfoil coefficients across a range of angle of attacks
# c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(xr, yr, alpha, re, iter=100, zeroinit=false, printdata=true, reinit=true)

# # a vs cl
# pltcl = plot(xlabel="Angle of Attack (deg)", ylabel="Lift Coefficient", legend=:best)
# read_coordinates("NASAcl.csv")
# plot!(x, y, label="Experimental Data", markershape=:square)
# plot!(alpha, c_l, label="XFoil Data")
# display(pltcl)

# # cl vs cd
# pltcd = plot(xlabel="Lift Coefficient", ylabel="Drag Coefficient", legend=:best)
# read_coordinates("NASAcd.csv")
# plot!(x, y, label="Experimental Data", markershape=:square)
# plot!(c_l, c_d, label="XFoil Data")
# display(pltcd)

# # cl vs cm
# pltcm = plot(xlabel="Lift Coefficient", ylabel="Moment Coefficient", legend=:best)
# read_coordinates("NASAcm.csv")
# plot!(x, y, label="Experimental Data", markershape=:square)
# plot!(c_l, c_d, label="XFoil Data")
# display(pltcm)

# ----------------------------------------------------------------------------

# function analyze_airfoil_re()
#     plt = plot(xlabel="Angle of Attack (deg)", ylabel="Lift Coefficient", legend=:best)
#     for i in [10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
#         re = i
#         # initialize outputs
#         n_a = length(alpha)
#         c_l = zeros(n_a)
#         c_d = zeros(n_a)
#         c_dp = zeros(n_a)
#         c_m = zeros(n_a)
#         converged = zeros(Bool, n_a)

#         # determine airfoil coefficients across a range of angle of attacks
#         for i = 1:n_a
#             c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
#         end

#         # plot results
#         if i == 1
#             plot!(plt, alpha, c_l, label="Re = $re")
#         else
#             plot!(plt, alpha, c_l, label="Re = $re")
#         end
#     end
#     return plt
# end


# plt, cl, cd, cm = analyze_airfoil_re()

# display(plt)

plt = plot(xlabel="Lift Coefficient", ylabel="Moment Coefficient", legend=:best)

# xX, yX = read_coordinates_X("naca2412.txt")
# load_coordinates(xX, yX)
# xr, yr = Xfoil.pane()
# for re in [1e5, 1e6, 1e7, 1e8, 1e9]
#     c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(xr, yr, alpha, re, iter=100, zeroinit=false, printdata=true, reinit=true)
#     plot!(alpha, c_l, label="Re = $re")
# end


function do_xfoil(filename, seriesname)
    xX, yX = read_coordinates_X(filename)
    load_coordinates(xX, yX)
    xr, yr = Xfoil.pane()
    c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(xr, yr, alpha, re; mach=mach, iter=100, zeroinit=true, printdata=true, reinit=false)
    plot!(alpha, c_m, label="XFoil $seriesname")
end

do_xfoil("naca2412.txt", "2412")

# do_xfoil("Thickness/naca2401.txt", "2401")
# do_xfoil("Thickness/naca2406.txt", "2406")
# do_xfoil("Thickness/naca2412.txt", "2412")
# do_xfoil("Thickness/naca2418.txt", "2418")
# do_xfoil("Thickness/naca2424.txt", "2424")
# do_xfoil("Thickness/naca2430.txt", "2430")

# do_xfoil("Max_Camber/naca0012.txt", "0012")
# do_xfoil("Max_Camber/naca1412.txt", "1412")
# do_xfoil("Max_Camber/naca2412.txt", "2412")
# do_xfoil("Max_Camber/naca3412.txt", "3412")
# do_xfoil("Max_Camber/naca4412.txt", "4412")
# do_xfoil("Max_Camber/naca5412.txt", "5412")
# do_xfoil("Max_Camber/naca7412.txt", "7412")
# do_xfoil("Max_Camber/naca9412.txt", "9412")

# do_xfoil("Camber_Position/naca2112.txt", "2112")
# do_xfoil("Camber_Position/naca2212.txt", "2212")
# do_xfoil("Camber_Position/naca2312.txt", "2312")
# do_xfoil("Camber_Position/naca2412.txt", "2412")
# do_xfoil("Camber_Position/naca2512.txt", "2512")
# do_xfoil("Camber_Position/naca2712.txt", "2712")
# do_xfoil("Camber_Position/naca2912.txt", "2912")

display(plt)

nothing
