#!/usr/bin/env julia
__precompile__(true)


module CellAdhesion

CellAdhesionFloat = Float32
CellAdhesionInt = Int32
export CellAdhesionFloat, CellAdhesionInt

using Plots

export Bond, Cluster
export BondModel, BondType, Slip, Catch
export SlipBondModel, CatchBondModel

include("bondmodels.jl")

include("definitions.jl")
include("utility.jl")
include("dynamics.jl")
include("processing.jl")


end # module
