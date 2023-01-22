
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

     junction = one_step(junction, model, convert(CellAdhesionFloat, 1.0))

     (typeof(junction) == Interface) && size(junction.history)==(12,2)

end

@test _check_one_step(tol)




function _check_junction_simulation(tol)


    model = Model(Dict("model"=>k_on_constant, "k_on_0"=>3e-3), 
                          Dict("model"=>k_off_slip, "k_off_0"=>3e-4, "f_1e"=>0.055), 
                          Dict("load"=>"local", "dt"=>0.01))
    junction = init_junction(50, model, history = false)

    stress_break_v = zeros(400)
    time_break_v = zeros(400)

    for sim = 1:1:400
        junction, stress_break_v[sim], time_break_v[sim], step = junction_simulation(junction, model, 0.08, max_steps = 500000)
    end
    stress_break_mean = mean(stress_break_v)
    time_break_mean = mean(time_break_v)


    model1 = Model(Dict("model"=>k_on_constant, "k_on_0"=>3e-3), 
                        Dict("model"=>k_off_slip, "k_off_0"=>3e-4, "f_1e"=>0.055), 
                        Dict("load"=>"global", "dt"=>0.01))
    junction1 = init_junction(50, model1, history = false)

    stress_break_v1 = zeros(400)
    time_break_v1 = zeros(400)

    for sim = 1:1:400
        junction1, stress_break_v1[sim], time_break_v1[sim], step = junction_simulation(junction1, model1, 0.08, max_steps = 500000)
    end
    stress_break_mean1 = mean(stress_break_v1)
    time_break_mean1 = mean(time_break_v1)

    (time_break_mean<time_break_mean1) && (time_break_mean<2.0) && (time_break_mean<10.0)

end

@test _check_junction_simulation(tol)