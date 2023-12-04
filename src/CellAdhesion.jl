#!/usr/bin/env julia
__precompile__(true)


module CellAdhesion

CellAdhesionFloat = Float32
CellAdhesionInt = Int32
export CellAdhesionFloat, CellAdhesionInt

using Plots

export Bond, Cluster, SlipBondModel


include("definitions.jl")
include("utility.jl")
include("dynamics.jl")
include("processing.jl")


end # module
