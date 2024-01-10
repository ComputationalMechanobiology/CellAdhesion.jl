using Test 
using Statistics
using Plots
using CellAdhesion



const tol = (eps(CellAdhesion.CellAdhesionFloat))^(0.125)

println("===============================================")
println("Testing CellAdhesion")
println("===============================================")
println("\n")


include("test_utility.jl")
include("test_dynamics.jl")
include("test_processing.jl")