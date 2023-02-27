export Bond, Interface, Model



mutable struct Bond
    state::Union{BitArray, Bool}      # False 0  = open, True 1 = closed
    k_on::CellAdhesionFloat
    k_off::CellAdhesionFloat
    f::CellAdhesionFloat
    history::Union{Array{CellAdhesionFloat}, Bool}
end



mutable struct Interface
    bonds::Vector{Bond}                               
    state::Bool                            # False = Open, True = closed
    const n::Integer
    const l::CellAdhesionFloat             # Distance between bonds
end


struct Model
    k_on::Dict
    k_off::Dict
    param::Dict
end




function Base.setproperty!(x::Interface, s::Symbol, new_x::Union{Vector{CellAdhesionFloat}, Array{CellAdhesionFloat}})

    for i = 1:1:x.n
      setfield!(x.bonds[i], s, new_x[i])
    end
  
  end

