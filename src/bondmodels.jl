abstract type BondType end

struct BondModel{T}
  p::NamedTuple   # Model parameters
  BondModel{T}(; kwargs...) where T = new{T}(map(CellAdhesionFloat, (; kwargs...)))
end

# Allow accessing parameters directly via dot syntax
Base.getproperty(b::BondModel, s::Symbol) = 
    getproperty(getfield(b, :p), s)

# Provide tab-completion and reflection of available properties
Base.propertynames(b::BondModel) = keys(getfield(b,:p))



# k_on should always be a constant. This could be specialised if needed, but what could it depend on?
function k_on(m::BondModel{T}) where T
  return m.p.k_on
end


# Each model needs a k_off

struct Slip <: BondType end
struct Catch <: BondType end


# No need to force the type of f, as this function is only called with this type in principle.
function k_off(m::BondModel{Slip}, f)
  return m.k_off_0 .* exp.(f ./ m.f_1e)
end

function k_off(m::BondModel{Catch}, f)
  return m.k_off_0s .* exp.(f ./ m.f_1es) + m.k_off_0c .* exp.(-f ./ m.f_1ec)
end


function SlipBondModel(k_on::NamedTuple, k_off::NamedTuple)
  return(BondModel{Slip}(k_on=k_on.k_on_0, k_off_0=k_off.k_off_0, f_0=k_off.f_0))
end
function CatchBondModel(k_on::NamedTuple, k_off::NamedTuple)
  return(BondModel{Catch}(k_on=k_on.k_on_0, k_off_0s=k_off.k_off_0s, f_1es=k_off.f_1es, k_off_0c=k_off.k_off_0c, f_1ec=k_off.f_1ec))
end


# # How to define the bond model

# m=BondModel{Slip}(k_on=1, k_off_0=2, f_0=3)


# # How to define a new bond type

# struct MyBond <: BondType end

# function k_off(m::BondModel{MyBond}, f)
#   return m.p.k_off_0 .* sinh.(f ./ m.p.f_1e)
# end
