
println("===============================================")
println("Testing processing.jl")
println("===============================================")

function _interface_bonds()

    model1 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.0), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict())
    model2 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict())
    model3 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.8), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), Dict())
    n = 3

    v1 = interface(n, 1.0, 15.0, true, model1)
    v2 = interface(n, 1.0, 10.0, true, model2)
    v3 = interface(n, 1.0, 15.0, true, model3)

    ((getproperty.(v1.u[:], :state) == zeros(n)) && (getproperty.(v2.u[:], :state) == ones(n)) && (typeof(getproperty.(v3.u[:], :state)) == BitVector) && (v1.f == 15.0) 
        && (v1.state == false) && (v2.state==true) && (v3.state == true))

end
@test _interface_bonds()



function _interface_cluster()

    model = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.2), Dict("model"=>k_off_slip, "k_off_0"=>0.8, "f_1e"=>1), Dict())
    v4 = interface([2, 3], [1.0, 0.1], 15.0, [true, true], model)

    (length(v4.u)==2) && (length(v4.u[1].u)==3) && (v4.f == 15.0)

end
@test _interface_cluster()






# function _check_one_step(tol)


#      model = model = Model(Dict("model"=>k_on_constant, "k_on_0"=>0.5), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1e4), Dict("load"=>"global", "dt"=>0.5))
#      junction = init_junction(12, model, history = true)

#      junction = one_step(junction, model, convert(CellAdhesionFloat, 1.0))

#      (typeof(junction) == Interface) && size(junction.history)==(12,2)

# end

# @test _check_one_step(tol)




# function _check_junction_simulation(tol)


#     model = Model(Dict("model"=>k_on_constant, "k_on_0"=>3e-3), 
#                           Dict("model"=>k_off_slip, "k_off_0"=>3e-4, "f_1e"=>0.055), 
#                           Dict("load"=>"local", "dt"=>0.01))
#     junction = init_junction(50, model, history = false)

#     stress_break_v = zeros(400)
#     time_break_v = zeros(400)

#     for sim = 1:1:400
#         junction, stress_break_v[sim], time_break_v[sim], step = junction_simulation(junction, model, 0.08, max_steps = 500000)
#     end
#     stress_break_mean = mean(stress_break_v)
#     time_break_mean = mean(time_break_v)


#     model1 = Model(Dict("model"=>k_on_constant, "k_on_0"=>3e-3), 
#                         Dict("model"=>k_off_slip, "k_off_0"=>3e-4, "f_1e"=>0.055), 
#                         Dict("load"=>"global", "dt"=>0.01))
#     junction1 = init_junction(50, model1, history = false)

#     stress_break_v1 = zeros(400)
#     time_break_v1 = zeros(400)

#     for sim = 1:1:400
#         junction1, stress_break_v1[sim], time_break_v1[sim], step = junction_simulation(junction1, model1, 0.08, max_steps = 500000)
#     end
#     stress_break_mean1 = mean(stress_break_v1)
#     time_break_mean1 = mean(time_break_v1)

#     (time_break_mean<time_break_mean1) && (time_break_mean<2.0) && (time_break_mean<10.0)

# end

# @test _check_junction_simulation(tol)