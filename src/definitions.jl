

"""
    Bond
    Bond(state::{Bool}, f::{CellAdhesionFloat}, model)

    `Bond` struct contains sate, force and model of a single bond

    # Fields

    - `state`: false (0)  = open, true (1) = closed
    - `f`: force
    - `t`: time
    - `model`: BondModel describing the bond binding-unbinding dynamics
"""
mutable struct Bond{T}
    state::Bool   # False 0  = open, True 1 = closed
    f::CellAdhesionFloat              # Force applied to the bond
    model::T                          # BondModel data type (we do not use pointers -> will see later if necessary) 
end

# 
"""
    Cluster(u::Vector{Bond,Cluster}, state::{Bool}, f::{CellAdhesionFloat}, f_model::{Symbol}, n::{CellAdhesionInt}, l::{CellAdhesionFloat})

"""
mutable struct Cluster{T}
    u::Vector{T}                     # Unit element       
    state::Bool                         # False 0 = Open, True 1 = closed
    f::CellAdhesionFloat                # Force applied to the Interface
    const f_model::Symbol
    const n::CellAdhesionInt
    const l::CellAdhesionFloat          # Distance between unit elements
end


