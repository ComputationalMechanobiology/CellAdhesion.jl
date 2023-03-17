
println("===============================================")
println("Testing processing.jl")
println("===============================================")

function _interface_bonds()

    model1 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.0), Dict("model"=>k_off_slip, "k_off_0"=>1.0, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    model2 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    model3 = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.8), Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    n = 3

    v1 = interface(n, 1.0, 15.0, false, model1)
    v2 = interface(n, 1.0, 10.0, false, model2)
    v3 = interface(n, 1.0, 15.0, false, model3)

    ((getproperty.(v1.u[:], :state) == zeros(n)) 
      && (getproperty.(v2.u[:], :state) == ones(n)) 
      && (typeof(getproperty.(v3.u[:], :state)) == BitVector) 
      && (v1.f == 15.0) 
      && (v1.state == false) 
      && (v2.state==true) 
      && (v3.state == true))

end
@test _interface_bonds()



function _interface_cluster()

    model = Model(Dict("model"=>force_global), Dict("model"=>k_on_constant, "k_on_0"=>0.2), Dict("model"=>k_off_slip, "k_off_0"=>0.8, "f_1e"=>1), Dict("dt"=>convert(CellAdhesionFloat,1e-2)))
    v4 = interface([2, 3], [1.0, 0.1], 15.0, [false, false], model)

    ((length(v4.u)==2) 
      && (length(v4.u[1].u)==3) 
      && (v4.f == 15.0))

end
@test _interface_cluster()



function _check_one_step(tol)

  model1 = Model(Dict("model"=>force_global), 
                Dict("model"=>k_on_constant, "k_on_0"=>1.0), 
                Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), 
                Dict("dt"=>convert(CellAdhesionFloat,1.0)))

  junction = interface(10, 1.0, 1.0, false, model1)
 
  one_step(junction, model1)

  (junction.state == true)  

end

@test _check_one_step(tol)




function _check_force_increment(tol)

  model1 = Model(Dict("model"=>force_global), 
                Dict("model"=>k_on_constant, "k_on_0"=>1.0), 
                Dict("model"=>k_off_slip, "k_off_0"=>0.2, "f_1e"=>1), 
                Dict("dt"=>convert(CellAdhesionFloat,1.0)))

  junction = interface(10, 1.0, 1.0, false, model1)

  force_increment(junction, model1, convert(CellAdhesionFloat,10.0))

  ((junction.f == 10.0) 
    && isapprox(sum(getfield.(junction.u, :f)), 10.0,atol=tol))

end

@test _check_force_increment(tol)




function _check_junction_simulation(tol)

  model1 = Model(Dict("model"=>force_global), 
                Dict("model"=>k_on_constant, "k_on_0"=>3e-3), 
                Dict("model"=>k_off_slip, "k_off_0"=>3e-4, "f_1e"=>0.055), 
                Dict("dt"=>convert(CellAdhesionFloat,0.01)))

  N = 20
  junction = interface(50, 1.0, 0.0, false, model1)
  stress_break_v = zeros(N)
  time_break_v = zeros(N)

  for sim = 1:1:N
    junction = interface(50, 1.0, 0.0, false, model1)
    stress_break_v[sim], time_break_v[sim], step = junction_simulation(junction, model1, 50*0.2, max_steps = 500000)
  end
  stress_break_mean = mean(stress_break_v)
  time_break_mean = mean(time_break_v)


  model2 = Model(Dict("model"=>force_local), 
                  Dict("model"=>k_on_constant, "k_on_0"=>3e-3), 
                  Dict("model"=>k_off_slip, "k_off_0"=>3e-4, "f_1e"=>0.055), 
                  Dict("dt"=>convert(CellAdhesionFloat,0.01)))
  
  

    stress_break_v1 = zeros(N)
    time_break_v1 = zeros(N)

    for sim = 1:1:N
      junction = interface(10, 1.0, 0.0, false, model2)
      stress_break_v1[sim], time_break_v1[sim], step = junction_simulation(junction, model2, 50*0.2, max_steps = 500000)
    end
    stress_break_mean1 = mean(stress_break_v1)
    time_break_mean1 = mean(time_break_v1)

    print(time_break_mean, "\n")
    print(time_break_mean1, "\n")

    (time_break_mean>time_break_mean1) #&& (time_break_mean<2.0) && (time_break_mean<10.0)

end

@test _check_junction_simulation(tol)