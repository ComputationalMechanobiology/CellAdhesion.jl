
println("===============================================")
println("Testing utility.jl")
println("===============================================")





function _update_state_bonds()

    model1 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.0), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict())
    model2 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict())
    model3 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.8), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), Dict())

    int_1 = interface(5, 1.0, 15.0, true, model1)
    int_2 = interface(5, 1.0, 15.0, true, model2)
    int_3 = interface(5, 1.0, 15.0, true, model3)

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)

    int_1.state == false && int_2.state == true && int_3.state == true

end
@test _update_state_bonds()


function _update_state_clusters()

    model1 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.0), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict())
    model2 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict())
    model3 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.8), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), Dict())


    int_1 = interface([2, 3], [1.0, 0.1], 15.0, [true, true], model1)
    int_2 = interface([2, 3], [1.0, 0.1], 15.0,[true, true], model2)
    int_3 = interface([2, 3], [1.0, 0.1], 15.0, [true, true], model3)

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)


    (int_1.state == false) && (int_2.state == true) && (int_3.state == true) && (int_1.u[1].state ==false) && (int_1.u[2].state ==false) && (int_2.u[1].state ==true)

end
@test _update_state_clusters()


