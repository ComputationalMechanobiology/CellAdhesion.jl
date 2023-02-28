export Bond, Interface, Model



mutable struct Bond
    state::Union{BitArray, Bool}      # False 0  = open, True 1 = closed
    k_on::CellAdhesionFloat
    k_off::CellAdhesionFloat
    f::CellAdhesionFloat
    history::Union{Array{CellAdhesionFloat}, Bool}
end



mutable struct Interface
    u::Union{Vector{Bond}, Vector{Interface}}            # Unit element (it can be a bond struct or an interface struct)       
    state::Bool                                          # False = Open, True = closed
    f::CellAdhesionFloat
    history::Union{Array{CellAdhesionFloat}, Bool} 
    const n::Integer
    const l::CellAdhesionFloat                           # Distance between bonds
end


struct Model
    k_on::Dict
    k_off::Dict
    param::Dict
end





