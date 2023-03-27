
println("===============================================")
println("Testing processing.jl")
println("===============================================")


function _check_KineticMonteCarlo_Cluster(tol)

  model1 = slip_model_init((k_on_0=0.2,), (k_off_0=1.0, f_1e=1))
  model2 = slip_model_init((k_on_0=0.8,), (k_off_0=0.0, f_1e=1))


  n = convert(CellAdhesionInt, 4)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 10.0)

  force_string = :force_global
  v1 = Cluster(Bond.([true,true,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), true, F, force_string, n, l)
  v2 = Cluster(Bond.([false,false,false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), false, F, force_string, n, l)

  update!(v1, convert(CellAdhesionFloat,10))
  update!(v2, convert(CellAdhesionFloat,10))


  (
    isapprox(getfield.(v1.u, :state), [false,false,false,false], atol=tol)
    && isapprox(getfield.(v2.u, :state), [true,true,true,true], atol=tol)
  )


end

@test _check_KineticMonteCarlo_Cluster(tol)


function _check_KineticMonteCarlo_Interface(tol)

  model1 = slip_model_init((k_on_0=0.2,), (k_off_0=1.0, f_1e=1))
  model2 = slip_model_init((k_on_0=0.8,), (k_off_0=0.0, f_1e=1))


  n = convert(CellAdhesionInt, 4)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 0.2)

  force_string = :force_global
  v1 = Cluster(Bond.([true,true,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([true,true,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_1 = Cluster([v1, v2], true, F, :force_global, convert(CellAdhesionInt, 2), l)
  setforce!(int_1, F)
  update!(int_1, convert(CellAdhesionFloat, 10))


  v3 = Cluster(Bond.([false,false,false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v4 = Cluster(Bond.([false,false,false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_2 = Cluster([v3, v4], true, F, :force_global, convert(CellAdhesionInt, 2), l)
  setforce!(int_2, F)
  update!(int_2, convert(CellAdhesionFloat, 10))


  (
    isapprox(getfield.(int_1.u[1].u, :state), [false,false,false,false], atol=tol)
    && isapprox(getfield.(int_1.u[2].u, :state), [false,false,false,false], atol=tol)
    && isapprox(getfield.(int_2.u[1].u, :state), [true,true,true, true], atol=tol)
    && isapprox(getfield.(int_2.u[2].u, :state), [true,true,true, true], atol=tol)
  )


end

@test _check_KineticMonteCarlo_Interface(tol)







# function _interface_bonds()

#   model1 = model_init((model=k_on_constant, k_on_0=0.0), (model=k_off_slip, k_off_0=1.0, f_1e=1))
#   model2 = model_init((model=k_on_constant, k_on_0=1.0), (model=k_off_slip, k_off_0=0.0, f_1e=1))
#   model3 = model_init((model=k_on_constant, k_on_0=0.8), (model=k_off_slip, k_off_0=0.2, f_1e=1))
    
#   n = 3
#   f_model = :force_global
#   v1 = interface(n, 1.0, model1, 15.0, f_model)
#   v2 = interface(n, 1.0, model2, 10.0, f_model)
#   v3 = interface(n, 1.0, model3, 15.0, f_model)

#   ((getproperty.(v1.u[:], :state) == zeros(n)) 
#     && (getproperty.(v2.u[:], :state) == ones(n)) 
#     && (typeof(getproperty.(v3.u[:], :state)) == BitVector) 
#     && (v1.f == 15.0)
#     && (v1.state == false) 
#     && (v2.state==true) 
#     && (v3.state == true))

# end
# @test _interface_bonds()



# function _interface_cluster()

#   model = model_init((model=k_on_constant, k_on_0=0.0), (model=k_off_slip, k_off_0=1.0, f_1e=1))

#   f_model1 = :force_global
#   f_model2 = :force_local
#   v = interface([2, 3], [1.0, 0.1], model, 15.0, [f_model1, f_model2])

#   ((length(v.u)==2) 
#     && (length(v.u[1].u)==3) 
#     && (v.f == 15.0))

# end
# @test _interface_cluster()



# function _check_one_step(tol)

#   model1 = model_init((model=k_on_constant, k_on_0=1.0), (model=k_off_slip, k_off_0=0.2, f_1e=1))

#   junction = interface(10, 1.0, model1, 1.0, :force_global)
 
#   one_step(junction, convert(CellAdhesionFloat, 1.0))

#   (junction.state == true)  

# end

# @test _check_one_step(tol)




# function _check_force_increment(tol)

#   model1 = model_init((model=k_on_constant, k_on_0=1.0), (model=k_off_slip, k_off_0=0.2, f_1e=1))

#   junction = interface(10, 1.0, model1, 1.0, :force_global)

#   force_increment(junction, convert(CellAdhesionFloat,10.0))

#   ((junction.f == 10.0) 
#     && isapprox(sum(getfield.(junction.u, :f)), 10.0,atol=tol))

# end

# @test _check_force_increment(tol)




# function _check_junction_simulation(tol)


#   model = model_init((model=k_on_constant, k_on_0=3e-3), (model=k_off_slip, k_off_0=3e-4, f_1e=0.055))

#   N = 20
#   #junction = interface(50, 1.0, model1, 0.0, :force_global)
#   stress_break_v = zeros(N)
#   time_break_v = zeros(N)

#   for sim = 1:1:N
#     junction = interface(50, 1.0, model, 0.0, :force_global)
#     stress_break_v[sim], time_break_v[sim], step = junction_simulation(junction, 50*0.2, 0.01, max_steps = 500000)
#   end
#   stress_break_mean = mean(stress_break_v)
#   time_break_mean = mean(time_break_v)


#   stress_break_v1 = zeros(N)
#   time_break_v1 = zeros(N)

#   for sim = 1:1:N
#     junction = interface(10, 1.0, model, 0.0, :force_local)
#     stress_break_v1[sim], time_break_v1[sim], step = junction_simulation(junction, 50*0.2, 0.01, max_steps = 500000)
#   end
#   stress_break_mean1 = mean(stress_break_v1)
#   time_break_mean1 = mean(time_break_v1)

#   # print(time_break_mean, "\n")
#   # print(time_break_mean1, "\n")

#   (time_break_mean>time_break_mean1) #&& (time_break_mean<2.0) && (time_break_mean<10.0)

# end

# @test _check_junction_simulation(tol)