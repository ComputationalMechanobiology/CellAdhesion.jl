
println("===============================================")
println("Testing processing.jl")
println("===============================================")




function _check_init_junction(tol)

    model = Model(Dict("model"=>k_on_constant, "k_on_0"=>0.1), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict("load"=>"global"))
    junction = init_junction(12, model)

    model1 = Model(Dict("model"=>k_on_constant, "k_on_0"=>0.1), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict("load"=>"global"))
    junction1 = init_junction(12, model1, history = true)

    (typeof(junction) == Interface)  && (junction1.history == junction1.v)


end
@test _check_init_junction(tol)




function _check_one_step(tol)


     model = model = Model(Dict("model"=>k_on_constant, "k_on_0"=>0.5), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1e4), Dict("load"=>"global", "dt"=>0.5))
     junction = init_junction(12, model, history = true)

     junction = one_step(junction, model, 1.0)

     (typeof(junction) == Interface) && size(junction.history)==(12,2)

end

@test _check_one_step(tol)