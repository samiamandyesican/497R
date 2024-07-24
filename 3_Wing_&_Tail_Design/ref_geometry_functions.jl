"""
    ref_area_trap(yle, chord)

Return the area of a wing given chord and leading edge y-coordinates.

# Arguments
- `yle::Vector`: 2-element vector of the root and tip y-coordinates for a wing leading edge.
- `chord::Vector`: 2-element vector of the root and tip chord-lengths. 

# Example
```jldoctest
julia> ref_area_trap([0.0, 1.0], [2.0, 1.6])
3.6
```
"""
function ref_area(yle::Vector, chord::Vector)
    area = 0.0

    for i in 2:length(yle)
        area += (chord[i] + chord[i-1]) * (yle[i] - yle[i-1]) # area of a trapezoid times 2 since the wing is symmetric
    end

    return area
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
    stepped_vec(min, max, iterations)

Return an n-element Vector{Int64} of equally spaced elements where n = iterations + 1 and min and max are the initial and final values respectively.

# Example
```jldoctest
julia> get_steps_num(3, 7, 8)
8-element Vector{Int64}:
 3
 4
 4
 5
 5
 6
 6
 7
```
"""
function stepped_vec(min, max, iterations::Int)

    # get step size
    step = (max - min) / (iterations-1)
    # initialize step matrix
    steps = zeros(iterations)
    steps[1] = min

    for i in 2:(iterations)
        steps[i] = min + (i - 1) * step
    end

    return steps

end