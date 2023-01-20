
println("===============================================")
println("Testing rates.jl")
println("===============================================")


function _check_k_off_slip(tol)

    junction = Interface(false, 4, [0,1,1,1], [0.1, 0.1, 0.1, 0.1], [0.5, 0.5, 0.5, 0.5], [1, 1, 1, 1], BitMatrix(undef,4,0))
    model = Model(Dict("model"=>"k_on_constant", "k_on_0"=>0.2), Dict("model"=>"k_off_slip", "k_off_0"=>0.1, "f_1e"=>1), Dict())
    junction = k_off_slip(junction, model)

    isapprox(junction.k_off, [0, 0.27182818, 0.27182818, 0.27182818], atol=tol) && (typeof(junction.k_off) == Vector{CellAdhesionFloat})

end
@test _check_k_off_slip(tol)


function _check_k_on_constant(tol)

    junction = Interface(false, 4, [1,0,0,1], [0.1, 0.1, 0.1, 0.1], [0.5, 0.5, 0.5, 0.5], [1, 1, 1, 1], BitMatrix(undef,4,0))
    model = Model(Dict("model"=>"k_on_constant", "k_on_0"=>0.2), Dict("model"=>"k_off_slip", "k_off_0"=>0.1, "f_1e"=>1), Dict())
    junction = k_on_constant(junction, model)

    isapprox(junction.k_on, [0,0.2,0.2,0], atol=tol) && (typeof(junction.k_on) == Vector{CellAdhesionFloat})

end
@test _check_k_on_constant(tol)

