"""
    run_VLM(chord, ns, nc, CDiff_run)

Given a chord distribution, number of spanwise panels for the wing (root to tip), 
number of chordwise panels, and whether to calculate far-field induced drag, 
return lift and induced drag coefficients and reference area. If CDiff = false, run final analysis and return in addition the 
surfaces, circulation with cooresponding y-coordinates, theoretical circulation and its corresponding y-coordinates, sd, min error, and max error

# Arguments
- `chord::Vector`: Vector of the chords at n+1 equally spaced points (in the y-direction) where the number of spanwise panels is n.
- `ns::Int`: Number of spanwise panels.
- `nc::Int`: Number of chord-wise panels.
- `CDiff_run::Bool`: Whether to calculate the coefficient of induced drag in the far-field

# Example
```jldoctest
julia> run_VLM([1.0, 1.0, 1.0], 2, 4, true)
(0.44201638298511275, 0.006320903931537388, 8.0)
julia> run_VLM([1.0, 1.0, 1.0], 2, 4, false)
(0.44201638298511275, 0.006320903931537388, 8.0, "CDiff was not calculated!", 8.0, Matrix{SurfacePanel{Float64}}[[SurfacePanel{Float64}(...0.01549525908845173, 0.08190029296198514))
```
"""
function run_VLM(chord::Vector{Float64}, ns::Int, nc::Int, CDiff_run::Bool)
    yle, zle, theta, phi, fc, bref = define_geometry(ns) # calculate geometry given ns for reference as to the number of elements in the geometry's vectors
    
    # objective
    ## given wing geometry, generate panels in preparation for VLM
    grid, surface = wing_to_surface_panels(-0.25*chord, yle, zle, chord, theta, phi, ns, nc; fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=mirror)
    surfaces = [surface]
    Sref = ref_area(yle, chord)
    cref = Sref / bref # caclulate reference chord
    ref = Reference(Sref, cref, bref, rref, Vinf) # create Reference struct for use by VLM
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric) # perform VLM steady_analysis 
    CF, CM = body_forces(system; frame=Wind()) # unpack 
    CD, CY, CL = CF # unpack

    if CDiff_run == true
        CDiff = far_field_drag(system) # calculate far-field drag in Trefftz plane
    else
        CDiff = "CDiff was not calculated!" # if CDiff is not run, generate error message and store it in CDiff
        circulation, circ_y_coor = extract_circulation(system, grid)
        max_circ = 2 * Vinf * Sref * CL / (bref * pi)  # formula for theoretical max circulation

        function get_ideal_circ(y)
            return max_circ * sqrt(1 - ( 2 * y / bref )^2) # formula for theoretical elliptic circulation distribution
        end

        ideal_circ_range = range(-bref/2, bref/2, 100) # generate points to calculate ideal circulation
        ideal_circ = get_ideal_circ.(ideal_circ_range) # calculate ideal circulation across range
        error = get_error(circulation, circ_y_coor, get_ideal_circ)

        return CL, CDiff, Sref, surfaces, circulation, circ_y_coor, ideal_circ, ideal_circ_range, error
    end

    return CL, CDiff, Sref
end


"""
    get_error(x::Vector, y::Vector, func)

Given a distribution x, its range y, and a function to compare against, return sd, min error, and max error

# Arguments
- `x::Vector`: Vector of the chords at n+1 equally spaced points (in the y-direction) where the number of spanwise panels is n.
- `y::Vector`: Number of spanwise panels.
- `func::Function`: Number of chord-wise panels.

# Example
```jldoctest
julia> function example(x)
           return x^2
       end
julia> get_error([1, 3, 2, 4], [1, 2, 3, 4], example)
(6.96419413859206, 1.0, 12.0)
```
"""
function get_error(x::Vector, y::Vector, func)
    # initialize variables for error functions
    n = length(x)
    total_variance = 0
    min_var = 0
    max_var = 0
    min_percent_error = 0
    max_percent_error = 0

    for i in 1:length(x)
        var = (x[i] - func(y[i]))^2 # variance
        total_variance += var
        if min_var == 0 # initialize minimum variance
            min_var = var   
        end 

        # find min and max variance
        if var > max_var
            max_var = var
        elseif var < min_var
            min_var = var
        end

    end

    for i in 1:length(x)
        perc_err = sqrt((x[i] - func(y[i]))^2) / func(y[i]) * 100 # variance
        if min_percent_error == 0 # initialize minimum variance
            min_percent_error = perc_err 
        end 

        # find min and max percent error
        if perc_err > max_percent_error
            max_percent_error = perc_err
        elseif perc_err < min_percent_error
            min_percent_error = perc_err
        end

    end

    sd = sqrt(total_variance / n) # standard deviation from variance
    min_error = sqrt(min_var) # min error from min variance
    max_error = sqrt(max_var) # max error from max variance
    min_perc_error = min_percent_error
    max_perc_error = max_percent_error
    error = sd, min_error, max_error, min_perc_error, max_perc_error
    return error
end


"""
    extract_circulation(system, grid)

Given the output 'system' from steady_analysis() and the corresponding yle, ns, and nc, return the circulation and the corresponding equally spaced y-coordinates.

# Arguments
- `system::System`: output of steady_analysis() from the VortexLattice.jl package.
- `grid::Array{Float64, 3}`: Ouput of wing_to_surface_panels() from the VortexLattice.jl package.

# Example
```jldoctest
julia> extract_circulation(system, grid)
[0.11105486781036324, 0.16203364037393395, 0.19539958142406819, 0.21977548214544987, 0.23731973434639814, 0.24980902196619076, 0.25771715208458146, 0.26163514978282976, 0.2616351497828298, 0.25771715208458135, 0.2498090219661906, 0.23731973434639808, 0.21977548214544973, 0.19539958142406819, 0.16203364037393395, 0.11105486781036336][-3.751908897254242, -3.254070851490458, -2.754325939316164, -2.2540164026121587, -1.7533439398577353, -1.2525016010391405, -0.7515333930966144, -0.2505233143829971, 0.2505233143829971, 0.7515333930966144, 1.2525016010391405, 1.7533439398577353, 2.2540164026121587, 2.754325939316164, 3.254070851490458, 3.751908897254242]
```
"""
function extract_circulation(system::System, grid)
    # combine all grid representations of surfaces into a single vector
    grids = [grid]

    # calculate lifting line geometry 
    ## r is the x, y, and z coordinates of each lifting line. 
    ## c is the chord length at each lifting line coordinate. 
    r, c = lifting_line_geometry(grids)

    # calculate lifting line coefficients
    cf, cm = lifting_line_coefficients(system, r, c; frame=Body())

    # unpack section lift coefficient
    cl = cf[1][end,:]

    function average_consecutives(x::Vector)
        n = length(x)-1
        # allocate new vector 
        newx = zeros(n)

        # find mean between consecutive elements of the vector and store in newx
        for i in 1:n
            newx[i] = (x[i+1]+x[i])/2
        end
        return newx
    end

    y_coor = average_consecutives(r[1][2,:]) # y coordinates spaced between each pair of consequtive lifting lines
    cvec = average_consecutives(c[1]) # chord at each y_coor

    circulation = @. cl * Vinf * cvec / 2 # formula for obtaining circulation of each lifting line from its c_l, Vinf, and chord length
    return circulation, y_coor
end


nothing