export Bond, Cluster, Model, SlipBondModel

mutable struct Bond{T}
    state::Bool   # False 0  = open, True 1 = closed
    f::CellAdhesionFloat              # Force applied to the bond
    model::T                          # BondModel data type (we do not use pointers -> will see later if necessary) 
end


mutable struct Cluster{T}
    u::Vector{T}                     # Unit element       
    state::Bool                         # False 0 = Open, True 1 = closed
    f::CellAdhesionFloat                # Force applied to the Interface
    const f_model::Symbol
    const n::CellAdhesionInt
    const l::CellAdhesionFloat          # Distance between unit elements
end


abstract type BondModel end

struct SlipBondModel <: BondModel
    k_on::NamedTuple
    k_off::NamedTuple

    function SlipBondModel(k_on::NamedTuple, k_off::NamedTuple)
        k_on_0 = convert(CellAdhesionFloat, k_on[:k_on_0])
        k_off_0 = convert(CellAdhesionFloat, k_off[:k_off_0])
        f_1e = convert(CellAdhesionFloat, k_off[:f_1e])
        new((k_on_0 = k_on_0,), (k_off_0 = k_off_0, f_1e = f_1e))
    end

end

# Interface = Cluster{Cluster}



# junction = Interface(...)

# setforce!(junction,f)
# update!(junction)

# setforce!(junction, f)
# # junction.f = 0
# # distribute_force(f_model, junction)
# # # for each element of u, setforce!(u[i], f/n)

# setforce!(bond,f)
# # bond.f=f

# update!(Bond) --> calculate k_on, k_off, set value of state (Monte-Carlo)

# update!(Cluster)
# # for each element of u, update!(u[i])
# # state = OR [u[i].state]





# read_state() of the cluster

