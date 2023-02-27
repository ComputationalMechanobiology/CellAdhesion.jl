using Test 
using Statistics
using CellAdhesion


const tol = (eps(CellAdhesion.CellAdhesionFloat))^(0.125)

println("===============================================")
println("Testing CellAdhesion")
println("===============================================")
println("|")
println("|")
println("|")

include("test_utility.jl")
#include("test_rates.jl")
#include("test_processing.jl")