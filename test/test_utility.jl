
println("===============================================")
println("Testing utility.jl")
println("===============================================")


function _check_slip_model_init()

    k_on_params = (k_on_0 = 0.2,)
    k_off_params = (k_off_0 = 0.8, f_1e = 1.0)

    model = slip_model_init(k_on_params, k_off_params)

    typeof(model) == SlipBondModel

end

@test _check_slip_model_init()

# function _update_state_bonds()

#     model1 = model_init((model=k_on_constant, k_on_0=0.0), (model=k_off_slip, k_off_0=1.0, f_1e=1))
#     model2 = model_init((model=k_on_constant, k_on_0=1.0), (model=k_off_slip, k_off_0=0.0, f_1e=1))
#     model3 = model_init((model=k_on_constant, k_on_0=0.8), (model=k_off_slip, k_off_0=0.2, f_1e=1))

#     int_1 = interface(5, 1.0, model1, 15.0, :force_global)
#     int_2 = interface(5, 1.0, model2, 15.0, :force_global)
#     int_3 = interface(5, 1.0, model3, 15.0, :force_global)

#     update_state(int_1)
#     update_state(int_2)
#     update_state(int_3)

#     ((int_1.state == false)
#       && (int_2.state == true)
#       && (int_3.state == true))

# end
# @test _update_state_bonds()


# function _update_state_clusters()

#     model1 = model_init((model=k_on_constant, k_on_0=0.0), (model=k_off_slip, k_off_0=1.0, f_1e=1))
#     model2 = model_init((model=k_on_constant, k_on_0=1.0), (model=k_off_slip, k_off_0=0.0, f_1e=1))
#     model3 = model_init((model=k_on_constant, k_on_0=0.8), (model=k_off_slip, k_off_0=0.2, f_1e=1))


#     int_1 = interface([2, 3], [1.0, 0.1], model1, 15.0, [:force_global, :force_global])
#     int_2 = interface([2, 3], [1.0, 0.1], model2, 15.0, [:force_global, :force_global])
#     int_3 = interface([2, 3], [1.0, 0.1], model3, 15.0, [:force_global, :force_global])

#     update_state(int_1)
#     update_state(int_2)
#     update_state(int_3)


#     ((int_1.state == false) 
#       && (int_2.state == true) 
#       && (int_3.state == true) 
#       && (int_1.u[1].state ==false) 
#       && (int_1.u[2].state ==false) 
#       && (int_2.u[1].state ==true))

# end
# @test _update_state_clusters()



# function _check_KineticMonteCarlo_Cluster(tol)

#   model1 = model_init((model=k_on_constant, k_on_0=0.2), (model=k_off_slip, k_off_0=1.0, f_1e=1))
#   model2 = model_init((model=k_on_constant, k_on_0=0.8), (model=k_off_slip, k_off_0=0.0, f_1e=1))

#   n = 5
#   l = 1
#   v1 = Cluster(Bond.([true,true,true,true,true], zeros(n), zeros(n), zeros(n), repeat([model1],n)), false, 10.0, :force_global, n, l)
#   update_state(v1)
#   force(v1)
#   v2 = Cluster(Bond.([false,false,false,false,false], zeros(n), zeros(n), zeros(n), repeat([model2],n)), false, 10.0, :force_global, n, l)
#   update_state(v2)
#   force(v2)
#   k_rate_junction(v2)

#   KineticMonteCarlo(v1, convert(CellAdhesionFloat,1e-2))
#   KineticMonteCarlo(v2, convert(CellAdhesionFloat,1e-2))


#   (
#     isapprox(getfield.(v1.u, :state), [false,false,false,false,false], atol=tol)
#     && isapprox(getfield.(v2.u, :state), [true,true,true,true,true], atol=tol)
#   )


# end

# _check_KineticMonteCarlo_Cluster(tol)



# function _check_KineticMonteCarlo_Interface(tol)

#   model1 = model_init((model=k_on_constant, k_on_0=0.2), (model=k_off_slip, k_off_0=1.0, f_1e=1))
#   model2 = model_init((model=k_on_constant, k_on_0=0.8), (model=k_off_slip, k_off_0=0.0, f_1e=1))
  
#   n = 5
#   l = 1
#   v1 = Cluster(Bond.([true,true,true, true, true], zeros(n), zeros(n), zeros(n), repeat([model1],n)), false, 0.0, :force_global, n, l)
#   v2 = Cluster(Bond.([true,true,true, true, true], zeros(n), zeros(n), zeros(n), repeat([model1],n)), false, 0.0, :force_global, n, l)
#   int_1 = Interface([v1, v2], false, 0.2, :force_global, 2, l)
#   update_state(int_1)
#   force(int_1)
#   k_rate_junction(int_1)


#   v3 = Cluster(Bond.([false,false,false, false, false], zeros(n), zeros(n), zeros(n), repeat([model2],n)), false, 0.0, :force_global, n, l)
#   v4 = Cluster(Bond.([false,false,false, false, false], zeros(n), zeros(n), zeros(n), repeat([model2],n)), false, 0.0, :force_global, n, l)
#   int_2 = Interface([v3, v4], false, 0.2, :force_global, 2, l)
#   update_state(int_2)
#   force(int_2)
#   k_rate_junction(int_2)


#   KineticMonteCarlo(int_1, convert(CellAdhesionFloat,1e-2))
#   KineticMonteCarlo(int_2, convert(CellAdhesionFloat,1e-2))


#   (
#     isapprox(getfield.(int_1.u[1].u, :state), [false,false,false,false,false], atol=tol)
#     && isapprox(getfield.(int_1.u[2].u, :state), [false,false,false,false,false], atol=tol)
#     && isapprox(getfield.(int_2.u[1].u, :state), [true,true,true, true, true], atol=tol)
#     && isapprox(getfield.(int_2.u[2].u, :state), [true,true,true, true, true], atol=tol)
#   )


# end

# _check_KineticMonteCarlo_Interface(tol)