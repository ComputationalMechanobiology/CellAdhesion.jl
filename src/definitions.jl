export Interface

"""
    Interface()

    `Interface` struct contains .



    # Fields

    - XX: 
"""

struct Interface
    state::Bool              # If 0 it is broken, if 1 it is still close
    n::Integer                  
    v::Vector{Integer}
    k_on::Vector{CellAdhesionFloat}
    k_off::Vector{CellAdhesionFloat}
    f::Vector{CellAdhesionFloat}
    param
end


