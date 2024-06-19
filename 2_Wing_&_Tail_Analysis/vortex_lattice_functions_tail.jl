localpath = @__DIR__
cd(localpath)
include("vortex_lattice_functions_wing_only.jl")


function tail_volume_h(yle, yle_h, chord, chord_h, translation)
    return ref_area_trap(yle_h, chord_h) * lever_h(yle, chord, yle_h, chord_h, translation) / (ref_area_trap(yle, chord) * ref_chord(chord))
end


function tail_volume_v(yle, zle_v, chord, chord_v, translation)
    return 0.5 * ref_area_trap(zle_v, chord_v) * lever_v(yle, chord, zle_v, chord_v, translation) / (ref_area_trap(yle, chord) * ref_span(yle))
end