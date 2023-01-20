export Interface, Model

"""
    Interface()

    `Interface` struct contains .



    # Fields

    - XX: 
"""

struct Interface
    state::Bool              # If 0 it is broken, if 1 it is still close
    n::Integer                  
    v::BitVector
    k_on::Vector{CellAdhesionFloat}
    k_off::Vector{CellAdhesionFloat}
    f::Vector{CellAdhesionFloat}
    history::BitMatrix
end


struct Model
    k_on::Dict
    k_off::Dict
    param::Dict
end

