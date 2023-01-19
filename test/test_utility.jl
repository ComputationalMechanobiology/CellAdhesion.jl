
println("===============================================")
println("Testing utility.jl")
println("===============================================")


function _init_bonds()
    n = 10

    v = init_bonds(n, convert(CellAdhesionFloat, 0.0))
    v == zeros(n)

    v = init_bonds(n,convert(CellAdhesionFloat, 1.0))
    v == ones(n)
end
@test _init_bonds()


function _check_state()

    junction1 = Interface(false, 5, [0,0,0,0,0], [], [], [], [])
    junction1 = check_state(junction1)

    junction2 = Interface(true, 5, [0,0,1,0,0], [], [], [], [])
    junction2 = check_state(junction2)

    junction1.state == true && junction2.state == false

end
@test _check_state()



function _check_link_k_off_slip(tol)

    junction = Interface(false, 4, [0,1,1,1], [0.1, 0.1, 0.1, 0.1], [0.5, 0.5, 0.5, 0.5], [1, 1, 1, 1], Dict("k_off_0"=>0.1, "f_1e"=>1))
    junction = link_k_off_slip(junction)

    isapprox(junction.k_off, [0, 0.27182818, 0.27182818, 0.27182818], atol=tol)

end
@test _check_link_k_off_slip(tol)



function _check_link_k_on_constant(tol)

    junction = Interface(false, 4, [1,0,0,1], [0.1, 0.1, 0.1, 0.1], [0.5, 0.5, 0.5, 0.5], [1, 1, 1, 1], Dict("k_off_0"=>0.1, "f_1e"=>1, "k_on_0"=>0.2))
    junction = link_k_on_constant(junction)

    isapprox(junction.k_on, [0,0.2,0.2,0], atol=tol)

end
@test _check_link_k_on_constant(tol)