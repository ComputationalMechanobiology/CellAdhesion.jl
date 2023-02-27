export Bond, Interface, Model



struct Bond
    state::Union{BitArray, Bool}      # False 0  = open, True 1 = closed
    k_on::CellAdhesionFloat
    k_off::CellAdhesionFloat
    f::CellAdhesionFloat
    history::Array{CellAdhesionFloat}
end



struct Interface
    bonds::Vector{Bond}             
    n::Integer                  
    state::Union{BitArray, Bool}    # False = Open, True = closed
    l::CellAdhesionFloat             # Distance between bonds
end


struct Model
    k_on::Dict
    k_off::Dict
    param::Dict
end

