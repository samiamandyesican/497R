#=
BYU FLOW Lab
Samuel Nasman
17 June 2024
functions for vortex lattice research
=#

"""
    ref_area_trap(yle, chord)

Return the area of a trapezoidal wing given chord and leading edge y-coordinate.

# Arguments
- `yle::Vector`: 2-element vector of the root and tip y-coordinates for a wing leading edge.
- `chord::Vector`: 2-element vector of the root and tip chord-lengths. 

# Example
```jldoctest
julia> ref_area_trap([0.0, 1.0], [2.0, 1.6])
3.6
```
"""
function ref_area_trap(yle, chord)
    return (chord[2] + chord[1]) * yle[2]
end


"""
    ref_chord(chord)

Return the reference chord given the root and tip chord-length as a two-element vector.

# Example
```jldoctest
julia> ref_chord([2.0, 1.6])
1.8
```
"""
function ref_chord(chord)
    return 0.5 * (chord[1] + chord[2])
end


"""
    ref_span(yle)

Return the reference span given the root and tip y-coordinates of the leading edge as a two-element vector.

# Example
```jldoctest
julia> ref_span([0.0, 7.5])
15.0
```
"""
function ref_span(yle)
    return 2 * yle[2]
end


"""
    center_of_g(chord)

Return the estimated center of gravity (quarter-reference-chord) given the root and tip chord-length as a two-element vector.

# Example
```jldoctest
julia> center_of_g([2.0, 1.6])
3-element Vector{Float64}:
 0.45
 0.0
 0.0
```
"""
function center_of_g(chord)
    return [0.25 * ref_chord(chord), 0.0, 0.0]
end


"""
    wing_efficiency(C_L, C_Di, yle, chord)

Return the inviscid span efficiency given the lift coefficient, drag coefficient, leading edge y-coordinates, and root and tip chord.

# Example
```jldoctest
julia> wing_efficiency(5.0, 2.0, [0.0. 7.5], [2.2, 1.8])
0.5305164769729845
```
"""
function wing_efficiency(C_L, C_Di, yle, chord)
    return C_L^2 / (π * aspect_ratio(yle, chord) * C_Di)
end


"""
    aspect_ratio(yle, chord)

Return the aspect ratio given the leading edge y-coordinates, and root and tip chord lengths as 2-element vectors.

# Example
```jldoctest
julia> aspect_ratio([0.0. 7.5], [2.2, 1.8])
7.5
```
"""
function aspect_ratio(yle, chord)
    return ref_span(yle)^2 / ref_area_trap(yle, chord)
end


"""
    force(CL, Sref, Vinf)

Return the lift given lift coefficient, wing reference area, and freestream velocity.

# Example
```jldoctest
julia> force(0.5, 20, 1)
5.0
```
"""
function force(CL, Sref, Vinf)
    return 0.5 * CL * Vinf^2 * Sref
end


"""
    get_step_float(min, max, iterations)

Return a list of equally spaced floats given a maximum value, minimum value, and number of iterations.

# Example
```jldoctest
julia> get_step_float(3.0, 4.0, 7)
7-element Vector{Float64}:
 3.0
 3.1666666666666665
 3.3333333333333335
 3.5
 3.6666666666666665
 3.833333333333333
 4.0
```
"""
function get_step_float(min, max, iterations)

    # get step size
    step = (max - min) / (iterations - 1)
    # initialize step matrix
    steps = zeros(iterations)
    steps[1] = min

    for i in 2:iterations
        steps[i] = min + (i - 1) * step
    end

    return steps

end

nothing