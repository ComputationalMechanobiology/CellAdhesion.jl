
println("===============================================")
println("Testing processing.jl")
println("===============================================")


function _check_KineticMonteCarlo_ClusterBonds(tol)

  model1 = SlipBondModel((k_on_0=0.2,), (k_off_0=1.0, f_1e=1))
  model2 = SlipBondModel((k_on_0=0.8,), (k_off_0=0.0, f_1e=1))


  n = convert(CellAdhesionInt, 4)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 10.0)

  force_string = :force_global
  v1 = Cluster(Bond.([true,true,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([false,false,false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v3 = Cluster(Bond.([false,false,false, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)


  setforce!(v1, F)
  update!(v1, convert(CellAdhesionFloat,10))
  setforce!(v1, F)

  setforce!(v2, F)
  update!(v2, convert(CellAdhesionFloat,10))
  setforce!(v2, F)

  setforce!(v3, F)
  update!(v3, convert(CellAdhesionFloat,10))
  setforce!(v3, F)


  (
    isapprox(getfield.(v1.u, :state), [false,false,false,false], atol=tol)
    && isapprox(getfield.(v2.u, :state), [true, true, true, true], atol=tol)
    && isapprox(getfield.(v3.u, :state), [true,true,true,true], atol=tol)
  )


end

@test _check_KineticMonteCarlo_ClusterBonds(tol)


function _check_KineticMonteCarlo_ClusterCluster(tol)

  model1 = SlipBondModel((k_on_0=0.2,), (k_off_0=1.0, f_1e=1))
  model2 = SlipBondModel((k_on_0=10.0,), (k_off_0=0.0, f_1e=1))


  n = convert(CellAdhesionInt, 4)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 10.0)

  force_string = :force_global
  v1 = Cluster(Bond.([true,true,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([true,true,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), :force_global, convert(CellAdhesionInt, 2), l)
  setforce!(int_1, F)
  update!(int_1, convert(CellAdhesionFloat, 1))
  setforce!(int_1, F)


  v3 = Cluster(Bond.([false,false,false, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v4 = Cluster(Bond.([false,false,false, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_2 = Cluster([v3, v4], true, convert(CellAdhesionFloat, 0.0), :force_global, convert(CellAdhesionInt, 2), l)
  setforce!(int_2, F)
  update!(int_2, convert(CellAdhesionFloat, 1))
  setforce!(int_2, F)


  v5 = Cluster(Bond.([false,false,false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v6 = Cluster(Bond.([false,false,false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_3 = Cluster([v5, v6], false, convert(CellAdhesionFloat, 0.0), :force_global, convert(CellAdhesionInt, 2), l)
  setforce!(int_3, F)
  update!(int_3, convert(CellAdhesionFloat, 1))
  setforce!(int_3, F)

  (
    isapprox(getfield.(int_1.u[1].u, :state), [false,false,false,false], atol=tol)
    && isapprox(getfield.(int_1.u[2].u, :state), [false,false,false,false], atol=tol)
    && isapprox(getfield.(int_2.u[1].u, :state), [true,true,true, true], atol=tol)
    && isapprox(getfield.(int_2.u[2].u, :state), [true,true,true, true], atol=tol)
    && isapprox(getfield.(int_3.u[1].u, :state), [true,true,true, true], atol=tol)
    && isapprox(getfield.(int_3.u[2].u, :state), [true,true,true, true], atol=tol)
  )


end

@test _check_KineticMonteCarlo_ClusterCluster(tol)




function _check_KineticMonteCarlo_Cluster2(tol)

  model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  n = convert(CellAdhesionInt, 4)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 60.0)

  force_string = :force_global
  v1 = Cluster(Bond.([true,false, true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([false,false, false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v3 = Cluster(Bond.([true,true, false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v4 = Cluster(Bond.([true,false, true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v5 = Cluster(Bond.([false,false, false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v6 = Cluster(Bond.([false,false, false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)

  c1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  c2 = Cluster([v3, v4], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  c3 = Cluster([v5, v6], false, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  
  int_1 = Cluster([c1, c2, c3], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 3), l)
  setforce!(int_1, F)
  update!(int_1, convert(CellAdhesionFloat, 1))
  setforce!(int_1, F)


  (
    isapprox(getfield.(int_1.u, :state), [true, true, true], atol=tol)
    && isapprox(getfield.(int_1.u[1].u, :state), [true, true], atol=tol)
    && isapprox(getfield.(int_1.u[2].u, :state), [true, true], atol=tol)
    && isapprox(getfield.(int_1.u[3].u, :state), [true, true], atol=tol)
    && isapprox(getfield.(int_1.u[1].u[1].u, :state), [true, true, true, true], atol=tol)
    && isapprox(getfield.(int_1.u[1].u[2].u, :state), [true,true,true, true], atol=tol)
    && isapprox(getfield.(int_1.u[2].u[1].u, :state), [true, true, true, true], atol=tol)
    && isapprox(getfield.(int_1.u[2].u[2].u, :state), [true, true, true, true], atol=tol)
    && isapprox(getfield.(int_1.u[3].u[1].u, :state), [true,true,true, true], atol=tol)
    && isapprox(getfield.(int_1.u[3].u[2].u, :state), [true,true,true, true], atol=tol)

  )


end

@test _check_KineticMonteCarlo_Cluster2(tol)



function _check_bond()

  model1 = SlipBondModel((k_on_0=0.0,), (k_off_0=1.0, f_1e=1))
  x1 = Bond(model1)
  model2 = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  x2 = Bond(model2)

  ((x1.state == false) && (x2.state == true) &&  typeof(x1) == Bond{SlipBondModel})

end
@test _check_bond()


function _check_cluster()

  n = convert(CellAdhesionInt, 5)
  l = convert(CellAdhesionFloat, 1.0)
  f_model1 = :force_global
  f_model2 = :force_local

  model1 = SlipBondModel((k_on_0=0.0,), (k_off_0=1.0, f_1e=1))
  x1 = Cluster(n, l, model1, f_model1)
  model2 = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  x2 = Cluster(n, l, model2, f_model2)

  (isapprox(getfield.(x1.u, :state), [false,false,false,false,false], atol=tol)
  && isapprox(getfield.(x2.u, :state), [true,true,true,true,true], atol=tol))

end

@test _check_cluster()


function _check_cluster()

  #n = convert(Vector{CellAdhesionInt}, [2,5])
  n = [2,5]
  l = [1.0, 0.1]
  f_model = [:force_global, :f_local]
  model1 = SlipBondModel((k_on_0=0.0,), (k_off_0=1.0, f_1e=1))
  x1 = Cluster(n, l, model1, f_model)

  n = [3, 2, 5]
  l = [1.0, 0.1, 0.01]
  f_model = [:force_global, :f_local, :f_global]
  model1 = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  x2 = Cluster(n, l, model1, f_model)

  (
    (x1.n == 2)
    && (typeof(x1.u[1].u[1]) == Bond{SlipBondModel})
    && (typeof(x1) == Cluster{Cluster})
    && (typeof(x2) == Cluster{Cluster})
    && (x2.n == 3)
    && (x2.u[1].n == 2)
    && (x2.u[1].u[1].u[1].state ==true)
    && (x2.state == true)
    
  )

end

@test _check_cluster()

function _check_runcluster(tol)

  model = SlipBondModel((k_on_0=3e-3,), (k_off_0=3e-4, f_1e=0.055))

  N = 20

  stress_break_v = zeros(N)
  time_break_v = zeros(N)

  for sim = 1:1:N
    x = Cluster(50, 1.0, model, :force_global)
    state, stress_break_v[sim], time_break_v[sim], step = runcluster(x, 50*0.2, 0.01, max_steps = 500000, state_check = true)
  end
  stress_break_mean = mean(stress_break_v)
  time_break_mean = mean(time_break_v)


  stress_break_v1 = zeros(N)
  time_break_v1 = zeros(N)

  for sim = 1:1:N
    x = Cluster(50, 1.0, model, :force_local)
    state, stress_break_v1[sim], time_break_v1[sim], step = runcluster(x, 50*0.2, 0.01, max_steps = 500000, state_check = true)
  end
  stress_break_mean1 = mean(stress_break_v1)
  time_break_mean1 = mean(time_break_v1)



  (time_break_mean>time_break_mean1) && (time_break_mean<20.0) && (time_break_mean>10.0)

end

@test _check_runcluster(tol)


function _check_runcluster_statecheck(tol)

  model = SlipBondModel((k_on_0=3e-3,), (k_off_0=3e-4, f_1e=0.055))

  N = 20

  stress_break_v = zeros(N)
  time_break_v = zeros(N)

  for sim = 1:1:N
    x = Cluster(50, 1.0, model, :force_global)
    state, stress_break_v[sim], time_break_v[sim], step = runcluster(x, 50*0.2, 0.01, max_steps = 2000, state_check = false)
  end
  stress_break_mean = mean(stress_break_v)
  time_break_mean = mean(time_break_v)


  stress_break_v1 = zeros(N)
  time_break_v1 = zeros(N)

  for sim = 1:1:N
    x = Cluster(50, 1.0, model, :force_local)
    state, stress_break_v1[sim], time_break_v1[sim], step = runcluster(x, 50*0.2, 0.01, max_steps = 1000, state_check = false)
  end
  stress_break_mean1 = mean(stress_break_v1)
  time_break_mean1 = mean(time_break_v1)

 (isapprox(time_break_mean, 20.0, atol=0.015)) && (isapprox(time_break_mean1, 10.0, atol=0.015))

end

@test _check_runcluster_statecheck(tol)