localpath = @__DIR__
cd(localpath)
include("vortex_lattice_functions_wing_only.jl")


"""
    lever_h(chord_h, translation, rref)

Return the length of the lever arm for the horizontal tail given horizontal tail chord, tail translation, and reference center.

# Arguments
- `translation::Vector`: 3-element vector of translation of the tail.
- `chord_h::Vector`: 2-element vector of the root and tip chord-lengths for the horizontal stabilizer. 
- `rref::Vector`: 3-element vector of reference center (generally the center of gravity).

# Example
```jldoctest
julia> lever_h([0.0, 1.0], [5.6, 0.0, 0.5], [0.7, 0.0, 0.0])
5.423098745182499
```
"""
function lever_h(chord_h, translation, rref)
    return sqrt((translation[1] + ref_chord(chord_h) - rref[1])^2 + (translation[2] - rref[2])^2 + (translation[3] - rref[3])^2)
end


"""
    mac_v(chord_v, xle_v, zle_v)

Return the coordinates for the mean aerodynamic center of the vertical stabilizer (prior to tail translation) given the geometry of the vertical stabilizer.

# Arguments
- `chord_v::Vector`: 2-element vector of the root and tip chord-lengths for the vertical stabilizer.
- `xle_v::Vector`: 2-element vector of the root and tip x-coordinates of the leading edge of the vertical stabiizer.
- `zle_v::Vector`: 2-element vector of the root and tip z-coordinates of the leading edge of the vertical stabiizer.

# Example
```jldoctest
julia> mac_v([1.0, 2.0], [0.0, 0.3], [0.2, 4.0])
3-element Vector{Float64}:
 0.9444444444444446
 0.0
 2.2222222222222223
```
"""
function mac_v(chord_v, xle_v, zle_v)
    # For the origin of these formulas Google "mean aerodynamic center" images and see the plot I generated (https://www.desmos.com/calculator/xprjrkuzzk)
    b = 2 * zle_v[2] # for this formula to work 'b' is twice the span of the vertical stabilizer.
    x = -chord_v[2] + (xle_v[2] + 2 * chord_v[2] + chord_v[1]) * ((chord_v[1] + 2 * chord_v[2]) / (3 * (chord_v[2] + chord_v[1])))
    z = b * (chord_v[1] + 2 * chord_v[2]) / (6 * (chord_v[2] + chord_v[1]))
    return [x, 0.0, z]
end


"""
    lever_v(chord_v, xle_v, zle_v, translation, rref)

Return the length of the lever arm for the vertical tail given vertical tail geometry, tail translation, and reference center.

# Arguments
- `chord_v::Vector`: 2-element vector of the root and tip chord-lengths for the vertical stabilizer. 
- `xle_v::Vector`: 2-element vector of the root and tip x-coordinates of the leading edge of the vertical stabiizer.
- `zle_v::Vector`: 2-element vector of the root and tip z-coordinates of the leading edge of the vertical stabiizer.
- `translation::Vector`: 3-element vector of translation of the tail.
- `rref::Vector`: 3-element vector of reference center (generally the center of gravity).

# Example
```jldoctest
julia> lever_v([1.0, 2.0], [0.0, 0.3], [0.2, 4.0], [5.6, 0.0, 0.5], [0.7, 0.0, 0.0])
6.447326941559425
```
"""
function lever_v(chord_v, xle_v, zle_v, translation, rref)
    mac = mac_v(chord_v, xle_v, zle_v)
    return sqrt((translation[1] + mac[1] - rref[1])^2 + (translation[2] - rref[2])^2 + (translation[3] + mac[3] - rref[3])^2)
end


"""
    tail_volume_h(yle, yle_h, chord, chord_h, translation, rref)

Return the horizontal tail volume given wing geometry, horizontal tail geometry, tail translation, and reference center.

# Arguments
- `yle::Vector`: 2-element vector of the root and tip y-coordinates of the leading edge of the wing.
- `yle_h::Vector`: 2-element vector of the root and tip y-coordinates of the leading edge of the horizontal stabiizer.
- `chord::Vector`: 2-element vector of the root and tip chord-lengths of the wing. 
- `chord_h::Vector`: 2-element vector of the root and tip chord-lengths for the horizontal stabilizer. 
- `translation::Vector`: 3-element vector of translation of the tail.
- `rref::Vector`: 3-element vector of reference center (generally the center of gravity).

# Example
```jldoctest
julia> tail_volume_h([0.0, 7.0], [0.0, 2.0], [2.5, 1.5], [2.0, 1.6], [5.6, 0.0, 0.5], [0.7, 0.0, 0.0])
0.8638239616580996
```
"""
function tail_volume_h(yle, yle_h, chord, chord_h, translation, rref)
    return ref_area_trap(yle_h, chord_h) * lever_h(chord_h, translation, rref) / (ref_area_trap(yle, chord) * ref_chord(chord))
end


"""
    tail_volume_v(yle, xle_v, zle_v, chord, chord_v, translation, rref)

Return the vertical tail volume given wing geometry, vertical tail geometry, tail translation, and reference center.

# Arguments
- `yle::Vector`: 2-element vector of the root and tip y-coordinates of the leading edge of the wing.
- `xle_v::Vector`: 2-element vector of the root and tip x-coordinates of the leading edge of the vertical stabiizer.
- `zle_v::Vector`: 2-element vector of the root and tip z-coordinates of the leading edge of the vertical stabiizer.
- `chord::Vector`: 2-element vector of the root and tip chord-lengths of the wing. 
- `chord_v::Vector`: 2-element vector of the root and tip chord-lengths for the vertical stabilizer. 
- `translation::Vector`: 3-element vector of translation of the tail.
- `rref::Vector`: 3-element vector of reference center (generally the center of gravity).

# Example
```jldoctest
julia> tail_volume_v([0.0, 7.0], [0.0, 0.8], [0.0, 4.0], [3.5, 2.5], [2.0, 1.6], [5.6, 0.0, 0.5], [0.7, 0.0, 0.0])
0.08139630689046873
```
"""
function tail_volume_v(yle, xle_v, zle_v, chord, chord_v, translation, rref)
    return 0.5 * ref_area_trap(zle_v, chord_v) * lever_v(chord_v, xle_v, zle_v, translation, rref) / (ref_area_trap(yle, chord) * ref_span(yle))
end


"""
    step_it(min_tran, max_tran, iterations)

Return step size (3-element vector) given minimum translation, maximum translation, and number of iterations

# Arguments
- `min_tran::Vector`: 3-element vector of the initial position for the tail.
- `max_tran::Vector`: 3-element vector of the final position for the tail.
- `iterations::Int`: Number of iterations (integer).

# Example
```jldoctest
julia> step_it([0.8, 0.0, 0.0], [10.0, 0.0, 3.0], 20)
3-element Vector{Float64}:
 0.4842105263157894
 0.0
 0.15789473684210525
```
"""
function step_it(min_tran, max_tran, iterations)
    return (max_tran - min_tran) ./ (iterations-1)
end


"""
    get_steps_vec(min, max, iterations)

Return 2xn matrix with n columns of equally spaced vectors where n = iterations and min and max are the initial and final vectors respectively.

# Arguments
- `min::Vector`: 2-element vector of the initial vector value.
- `max::Vector`: 2-element vector of the final vector value.
- `iterations::Int`: Number of iterations (integer).

# Example
```jldoctest
julia> get_steps_vec([1.0, 2.0], [10.0, 21.0], 5)
2×5 Matrix{Float64}:
 1.0  3.25   5.5   7.75  10.0
 2.0  6.75  11.5  16.25  21.0
```
"""
function get_steps_vec(min, max, iterations)

    # get step size
    step = (max - min) / (iterations - 1)
    # initialize step matrix
    steps = zeros(2, iterations)
    steps[:, 1] = min

    for i in 2:iterations
        steps[:, i] = min + (i - 1) * step
    end

    return steps

end


"""
    get_steps_num(min, max, iterations)

Return an n-element Vector{Int64} of equally spaced elements rounded to integers where n = iterations and min and max are the initial and final values respectively.

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
function get_steps_num(min, max, iterations)

    # get step size
    step = (max - min) / (iterations - 1)
    # initialize step matrix
    steps = zeros(Int, iterations)
    steps[1] = round(Int, min)

    for i in 2:iterations
        steps[i] = round(Int, min + (i - 1) * step)
    end

    return steps

end


nothing