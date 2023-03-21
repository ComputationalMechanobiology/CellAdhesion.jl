module CellAdhesion

CellAdhesionFloat = Float32
CellAdhesionInt = Int32
export CellAdhesionFloat, CellAdhesionInt

using Plots



include("definitions.jl")
include("utility.jl")
include("dynamics.jl")
include("processing.jl")

end # module
