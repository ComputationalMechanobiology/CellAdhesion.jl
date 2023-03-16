module CellAdhesion

CellAdhesionFloat = Float32
export CellAdhesionFloat

using Plots



include("definitions.jl")
include("utility.jl")
include("dynamics.jl")
include("processing.jl")

end # module
