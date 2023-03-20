
println("===============================================")
println("Testing utility.jl")
println("===============================================")


function _update_state_bonds()

    model1 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.0), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    model2 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    model3 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.8), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))

    int_1 = interface(5, 1.0, 15.0, false, model1)
    int_2 = interface(5, 1.0, 15.0, false, model2)
    int_3 = interface(5, 1.0, 15.0, false, model3)

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)

    ((int_1.state == false)
      && (int_2.state == true)
      && (int_3.state == true))

end
@test _update_state_bonds()


function _update_state_clusters()

    model1 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.0), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    model2 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    model3 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.8), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))


    int_1 = interface([2, 3], [1.0, 0.1], 15.0, [false, false], model1)
    int_2 = interface([2, 3], [1.0, 0.1], 15.0,[false, false], model2)
    int_3 = interface([2, 3], [1.0, 0.1], 15.0, [false, false], model3)

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)


    ((int_1.state == false) 
      && (int_2.state == true) 
      && (int_3.state == true) 
      && (int_1.u[1].state ==false) 
      && (int_1.u[2].state ==false) 
      && (int_2.u[1].state ==true))

end
@test _update_state_clusters()



function _check_KineticMonteCarlo_unit(tol)

  model1 = Model(Dict("model"=>force_global), 
                  Dict("model"=>k_on_constant, "k_on_0"=>0.2), 
                  Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), 
                  Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
  model2 = Model(Dict("model"=>force_global), 
                  Dict("model"=>k_on_constant, "k_on_0"=>0.8), 
                  Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), 
                  Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
  n = 5
  l = 1
  v1 = Interface(Bond.([true,true,true,true,true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)
  update_state(v1)
  force(v1, model1)
  v2 = Interface(Bond.([false,false,false,false,false], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)
  update_state(v2)
  force(v2, model1)
  k_rate_junction(v2, model1)


  KineticMonteCarlo_unit(v1, model1.param["dt"])
  KineticMonteCarlo_unit(v2, model1.param["dt"])


  (
    isapprox(getfield.(v1.u, :state), [false,false,false,false,false], atol=tol)
    && isapprox(getfield.(v2.u, :state), [true,true,true,true,true], atol=tol)
  )


end

_check_KineticMonteCarlo_unit(tol)



function _check_KineticMonteCarlo(tol)

  model1 = Model(Dict("model"=>force_global), 
                  Dict("model"=>k_on_constant, "k_on_0"=>0.2), 
                  Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), 
                  Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
  model2 = Model(Dict("model"=>force_global), 
                  Dict("model"=>k_on_constant, "k_on_0"=>0.8), 
                  Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), 
                  Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
  
  n = 5
  l = 1
  v1 = Interface(Bond.([true,true,true, true, true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
  v2 = Interface(Bond.([true,true,true, true, true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
  int_1 = Interface([v1, v2], false, 0.2, false, 2, l)
  update_state(int_1)
  force(int_1, model1)
  k_rate_junction(int_1, model1)


  v3 = Interface(Bond.([false,false,false, false, false], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
  v4 = Interface(Bond.([false,false,false, false, false], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
  int_2 = Interface([v3, v4], false, 0.2, false, 2, l)
  update_state(int_2)
  force(int_2, model2)
  k_rate_junction(int_2, model2)


  KineticMonteCarlo(int_1, model1)
  KineticMonteCarlo(int_2, model2)


  (
    isapprox(getfield.(int_1.u[1].u, :state), [false,false,false,false,false], atol=tol)
    && isapprox(getfield.(int_1.u[2].u, :state), [false,false,false,false,false], atol=tol)
    && isapprox(getfield.(int_2.u[1].u, :state), [true,true,true, true, true], atol=tol)
    && isapprox(getfield.(int_2.u[2].u, :state), [true,true,true, true, true], atol=tol)
  )


end

_check_KineticMonteCarlo(tol)