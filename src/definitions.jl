export Bond, Interface, Model


struct Model
    f::Dict
    k_on::Dict
    k_off::Dict
    param::Dict

end


mutable struct Bond
    state::Union{BitArray, Bool}      # False 0  = open, True 1 = closed
    k_on::CellAdhesionFloat
    k_off::CellAdhesionFloat
    f::CellAdhesionFloat              # Force applied to the bond
    model::Base.RefValue{Model}       # Pointer to model 
end



mutable struct Interface
    u::Union{Vector{Bond}, Vector{Interface}}            # Unit element (it can be a bond struct or an interface struct)       
    state::Bool                                          # False = Open, True = closed
    f::CellAdhesionFloat                                 # Force applied to the Interface
    const n::CellAdhesionInt
    const l::CellAdhesionFloat                           # Distance between bonds
end





# set_force!() change force on mutable interface

# read_state!() of the cluster

