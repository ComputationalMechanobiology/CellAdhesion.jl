
using JSON

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

  ((x1.state == false) && (x2.state == true) &&  typeof(x1) == Bond{BondModel{Slip}})

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


function _check_cluster_v2()

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
    && (typeof(x1.u[1].u[1]) == Bond{BondModel{Slip}})
    && (typeof(x1) == Cluster{Cluster})
    && (typeof(x2) == Cluster{Cluster})
    && (x2.n == 3)
    && (x2.u[1].n == 2)
    && (x2.u[1].u[1].u[1].state ==true)
    && (x2.state == true)
    
  )

end

@test _check_cluster_v2()

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


function _check_runcluster_json_output()

  model = SlipBondModel((k_on_0=0.0,), (k_off_0=0.0, f_1e=1.0))

  bonds = Bond.(
    [true, true],
    convert(Vector{CellAdhesionFloat}, zeros(2)),
    repeat([model], 2),
  )

  cluster = Cluster(
    bonds,
    true,
    convert(CellAdhesionFloat, 0.0),
    :force_global,
    convert(CellAdhesionInt, 2),
    convert(CellAdhesionFloat, 1.0),
  )

  force_history = [1.0, 2.0, 3.0]
  json_path = "test_runcluster_states.json"

  state, break_force, break_time, steps = runcluster(cluster, force_history, 0.5, json_path; max_steps = 3, verbose = false)

  saved = JSON.parsefile(json_path)

  cluster_ok = (
    saved["cluster"]["type"] == "cluster"
    && saved["cluster"]["f_model"] == "force_global"
    && saved["cluster"]["n"] == 2
    && isapprox(saved["cluster"]["l"], 1.0)
    && length(saved["cluster"]["children"]) == 2
  )

  force_ok = all(isapprox.(saved["force"], [1.0, 2.0, 3.0]))
  dt_ok = isapprox(saved["dt"], 0.5)

  time_key_map = Dict(
    parse(Float64, k) => k
    for k in keys(saved)
    if k != "cluster" && k != "force" && k != "dt"
  )

  times_ok = (
    length(time_key_map) == 3
    && haskey(time_key_map, 0.0)
    && haskey(time_key_map, 0.5)
    && haskey(time_key_map, 1.0)
  )

  t0_force_tree = saved[time_key_map[0.0]]
  t05_force_tree = saved[time_key_map[0.5]]
  t1_force_tree = saved[time_key_map[1.0]]

  force_trees_ok = (
    isapprox(t0_force_tree[1], 1.0)
    && all(isapprox.(t0_force_tree[2], [0.5, 0.5]))
    && isapprox(t05_force_tree[1], 2.0)
    && all(isapprox.(t05_force_tree[2], [1.0, 1.0]))
    && isapprox(t1_force_tree[1], 3.0)
    && all(isapprox.(t1_force_tree[2], [1.5, 1.5]))
  )

  run_ok = (state == true) && isapprox(break_force, convert(CellAdhesionFloat, 3.0)) && isapprox(break_time, convert(CellAdhesionFloat, 1.5)) && (steps == 3)

  return run_ok && cluster_ok && force_ok && dt_ok && times_ok && force_trees_ok

end

@test _check_runcluster_json_output()


function _check_bond_state_force()

    model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
    force_string = :force_global

    closed_bond = Bond(true, convert(CellAdhesionFloat, 2.5), model)
    open_bond = Bond(false, convert(CellAdhesionFloat, 7.5), model)

    bond_states_closed, bond_forces_closed = bond_state_force(closed_bond)
    bond_states_open, bond_forces_open = bond_state_force(open_bond)

    flat_cluster = Cluster(
        Bond.([true, false, true], convert(Vector{CellAdhesionFloat}, [1.0, 2.0, 3.0]), repeat([model], 3)),
        true,
        convert(CellAdhesionFloat, 0.0),
        force_string,
        convert(CellAdhesionInt, 3),
        convert(CellAdhesionFloat, 1.0),
    )

    nested_cluster = Cluster(
        [
            flat_cluster,
            Cluster(
                Bond.([false, true], convert(Vector{CellAdhesionFloat}, [4.0, 5.0]), repeat([model], 2)),
                true,
                convert(CellAdhesionFloat, 0.0),
                force_string,
                convert(CellAdhesionInt, 2),
                convert(CellAdhesionFloat, 1.0),
            ),
        ],
        true,
        convert(CellAdhesionFloat, 0.0),
        force_string,
        convert(CellAdhesionInt, 2),
        convert(CellAdhesionFloat, 1.0),
    )

    flat_states, flat_forces = bond_state_force(flat_cluster)
    nested_states, nested_forces = bond_state_force(nested_cluster)

    expected_flat_states = Bool[true, false, true]
    expected_flat_forces = convert(Vector{CellAdhesionFloat}, [1.0, NaN, 3.0])
    expected_nested_states = Bool[true, false, true, false, true]
    expected_nested_forces = convert(Vector{CellAdhesionFloat}, [1.0, NaN, 3.0, NaN, 5.0])

    ((bond_states_closed == Bool[true])
     && (bond_forces_closed == convert(Vector{CellAdhesionFloat}, [2.5]))
     && (bond_states_open == Bool[false])
     && (length(bond_forces_open) == 1)
     && isnan(bond_forces_open[1])
     && (flat_states == expected_flat_states)
     && all(isequal.(flat_forces, expected_flat_forces))
     && (nested_states == expected_nested_states)
     && all(isequal.(nested_forces, expected_nested_forces)))

end

@test _check_bond_state_force()

function _check_bond_state_force_nested()
    model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
    force_string = :force_global

    flat_cluster = Cluster(
        Bond.([true, false, true], convert(Vector{CellAdhesionFloat}, [1.0, 2.0, 3.0]), repeat([model], 3)),
        true,
        convert(CellAdhesionFloat, 0.0),
        force_string,
        convert(CellAdhesionInt, 3),
        convert(CellAdhesionFloat, 1.0),
    )

    nested_cluster = Cluster(
        [
            flat_cluster,
            Cluster(
                Bond.([false, true], convert(Vector{CellAdhesionFloat}, [4.0, 5.0]), repeat([model], 2)),
                true,
                convert(CellAdhesionFloat, 0.0),
                force_string,
                convert(CellAdhesionInt, 2),
                convert(CellAdhesionFloat, 1.0),
            ),
        ],
        true,
        convert(CellAdhesionFloat, 0.0),
        force_string,
        convert(CellAdhesionInt, 2),
        convert(CellAdhesionFloat, 1.0),
    )

    flat_nested = bond_state_force(flat_cluster; output = :nested, time = 0.5)
    nested_nested = bond_state_force(nested_cluster; output = :nested, time = 1.3)

    @test haskey(flat_nested, 0.5)
    flat_force_tree = flat_nested[0.5]
    @test isapprox(flat_force_tree[1], convert(CellAdhesionFloat, 0.0))
    @test all(isequal.(flat_force_tree[2], convert(Vector{CellAdhesionFloat}, [1.0, NaN, 3.0])))

    @test haskey(nested_nested, 1.3)
    nested_force_tree = nested_nested[1.3]
    @test isapprox(nested_force_tree[1], convert(CellAdhesionFloat, 0.0))
    child_forces = nested_force_tree[2]
    @test isapprox(child_forces[1][1], convert(CellAdhesionFloat, 0.0))
    @test all(isequal.(child_forces[1][2], convert(Vector{CellAdhesionFloat}, [1.0, NaN, 3.0])))
    @test isapprox(child_forces[2][1], convert(CellAdhesionFloat, 0.0))
    @test all(isequal.(child_forces[2][2], convert(Vector{CellAdhesionFloat}, [NaN, 5.0])))
    return true
end

@test _check_bond_state_force_nested()


function _check_save_cluster_state()
    model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
    c = Cluster(3, 1.0, model, :force_global)
    c.u[1].state = true; c.u[1].f = 2.0
    c.u[2].state = false; c.u[2].f = 2.0
    c.u[3].state = true; c.u[3].f = 3.0
    c.state = true
    c.f = convert(CellAdhesionFloat, 5.0)

    state_dict1 = bond_state_force(c; output = :nested, time = 0.0)
    state_dict2 = bond_state_force(c; output = :nested, time = 1.0)
    state_dict3 = bond_state_force(c; output = :nested, time = 2.0)
    states = [state_dict1, state_dict2, state_dict3]

    force_history = [1.0, 2.0, 3.0]
    dt = 1.0
    tmp1 = "test_state.json"
    save_cluster_state(c, force_history, states, tmp1, dt)

    saved = JSON.parsefile(tmp1; allownan = true, nan = "nan")

    cluster_ok = (
      saved["cluster"]["type"] == "cluster"
      && saved["cluster"]["f_model"] == "force_global"
      && saved["cluster"]["n"] == 3
      && isapprox(saved["cluster"]["l"], 1.0)
      && length(saved["cluster"]["children"]) == 3
    )

    force_ok = all(isapprox.(saved["force"], force_history))
    dt_ok = isapprox(saved["dt"], dt)

    time_key_map = Dict(
      parse(Float64, k) => k
      for k in keys(saved)
      if k != "cluster" && k != "force" && k != "dt"
    )

    times_ok = (
      length(time_key_map) == 3
      && haskey(time_key_map, 0.0)
      && haskey(time_key_map, 1.0)
      && haskey(time_key_map, 2.0)
    )

    t0_force_tree = saved[time_key_map[0.0]]
    t1_force_tree = saved[time_key_map[1.0]]
    t2_force_tree = saved[time_key_map[2.0]]

    trees_ok = (
      isapprox(t0_force_tree[1], 5.0)
      && isapprox(t0_force_tree[2][1], 2.0)
      && isnan(t0_force_tree[2][2])
      && isapprox(t0_force_tree[2][3], 3.0)
      && isapprox(t1_force_tree[1], 5.0)
      && isapprox(t1_force_tree[2][1], 2.0)
      && isnan(t1_force_tree[2][2])
      && isapprox(t1_force_tree[2][3], 3.0)
      && isapprox(t2_force_tree[1], 5.0)
      && isapprox(t2_force_tree[2][1], 2.0)
      && isnan(t2_force_tree[2][2])
      && isapprox(t2_force_tree[2][3], 3.0)
    )

    return cluster_ok && force_ok && dt_ok && times_ok && trees_ok
end

@test _check_save_cluster_state()


function _check_save_load_cluster_parameters()
  model = SlipBondModel((k_on_0=1.1,), (k_off_0=0.2, f_1e=3.3))

  child_a = Cluster(2, 0.5, model, :force_global)
  child_b = Cluster(2, 0.5, model, :force_local)

  root_cluster = Cluster(
    [child_a, child_b],
    true,
    convert(CellAdhesionFloat, 0.0),
    :force_global,
    convert(CellAdhesionInt, 2),
    convert(CellAdhesionFloat, 1.5),
  )

  tmp = "test_cluster_params.json"
  save_cluster_parameters(root_cluster, tmp)
  loaded = load_cluster_parameters(tmp)

  root_ok = (
    loaded.n == 2
    && loaded.f_model == :force_global
    && isapprox(loaded.l, convert(CellAdhesionFloat, 1.5))
  )

  children_ok = (
    length(loaded.u) == 2
    && loaded.u[1] isa Cluster
    && loaded.u[2] isa Cluster
    && loaded.u[1].n == 2
    && loaded.u[2].n == 2
    && loaded.u[1].f_model == :force_global
    && loaded.u[2].f_model == :force_local
  )

  bonds_ok = (
    all(b -> b isa Bond{BondModel{Slip}}, loaded.u[1].u)
    && all(b -> b isa Bond{BondModel{Slip}}, loaded.u[2].u)
  )

  p1 = getfield(loaded.u[1].u[1].model, :p)
  p2 = getfield(loaded.u[2].u[2].model, :p)
  params_ok = (
    isapprox(p1.k_on, convert(CellAdhesionFloat, 1.1))
    && isapprox(p1.k_off_0, convert(CellAdhesionFloat, 0.2))
    && isapprox(p1.f_1e, convert(CellAdhesionFloat, 3.3))
    && isapprox(p2.k_on, convert(CellAdhesionFloat, 1.1))
    && isapprox(p2.k_off_0, convert(CellAdhesionFloat, 0.2))
    && isapprox(p2.f_1e, convert(CellAdhesionFloat, 3.3))
  )

  reset_ok = (
    all(b -> (b.state == false) && isnan(b.f), loaded.u[1].u)
    && all(b -> (b.state == false) && isnan(b.f), loaded.u[2].u)
  )

  return root_ok && children_ok && bonds_ok && params_ok && reset_ok
end

@test _check_save_load_cluster_parameters()
