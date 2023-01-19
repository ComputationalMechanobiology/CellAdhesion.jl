using Revise
using Test 
using CellAdhesion

const tol = (eps(CellAdhesion.CellAdhesionFloat))^(0.125)

println("===============================================")
println("Testing CellAdhesion")
println("===============================================")
println("|")


include("test_utility.jl")