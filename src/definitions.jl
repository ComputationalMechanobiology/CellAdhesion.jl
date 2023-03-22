export Bond, Cluster, Interface, Model


struct Model
    k_on::NamedTuple
    k_off::NamedTuple
end


mutable struct Bond
    state::Union{BitArray, Bool}      # False 0  = open, True 1 = closed
    k_on::CellAdhesionFloat
    k_off::CellAdhesionFloat
    f::CellAdhesionFloat              # Force applied to the bond
    model::Base.RefValue{Model}       # Pointer to model 
end



mutable struct Cluster
    u::Vector{Bond}                     # Unit element       
    state::Bool                         # False = Open, True = closed
    f::CellAdhesionFloat                # Force applied to the Interface
    const f_model::Symbol
    const n::CellAdhesionInt
    const l::CellAdhesionFloat          # Distance between bonds
end


mutable struct Interface
    u::Vector{Cluster}                     # Unit element       
    state::Bool                         # False = Open, True = closed
    f::CellAdhesionFloat                # Force applied to the Interface
    const f_model::Symbol
    const n::CellAdhesionInt
    const l::CellAdhesionFloat          # Distance between bonds
end



# set_force!() change force on mutable interface

# read_state!() of the cluster

