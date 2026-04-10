
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


function _check_save_cluster_to_json()

  model = SlipBondModel((k_on_0=1.0,), (k_off_0=2.0, f_1e=3.0))
  cluster = Cluster(5, 1.0, model, :force_global)
  save_cluster_to_json(cluster, "cluster_params.json")

  loaded_cluster = load_from_json("cluster_params.json", false)

  (
    loaded_cluster isa Cluster
    && loaded_cluster.n == 5
    && isapprox(loaded_cluster.l, 1.0)
    && loaded_cluster.f_model == :force_global
    && all(b -> b isa Bond, loaded_cluster.u)
    && all(isapprox.([b.state for b in loaded_cluster.u], falses(5), atol=tol))
    && all(isnan.([b.f for b in loaded_cluster.u]))
    && all(b -> b.model isa BondModel{Slip}, loaded_cluster.u)
    && all(isapprox.([b.model.k_on for b in loaded_cluster.u], ones(5), atol=tol))
    && all(isapprox.([b.model.k_off_0 for b in loaded_cluster.u], fill(2.0, 5), atol=tol))
    && all(isapprox.([b.model.f_1e for b in loaded_cluster.u], fill(3.0, 5), atol=tol))
  )

end

@test _check_save_cluster_to_json()


function _check_load_from_json()

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

  # test case where force is not loaded
  cluster = load_from_json(json_path, false)
  no_force_ok = (
    cluster isa Cluster
    && cluster.n == 2
    && isapprox(cluster.l, 1.0)
    && cluster.f_model == :force_global
    && all(b -> b isa Bond, cluster.u)
    && all(isapprox.([b.state for b in cluster.u], falses(2), atol=tol))
    && all(isnan.([b.f for b in cluster.u]))
    && all(b -> b.model isa BondModel{Slip}, cluster.u)
    && all(isapprox.([b.model.k_on for b in cluster.u], zeros(2), atol=tol))
    && all(isapprox.([b.model.k_off_0 for b in cluster.u], zeros(2), atol=tol))
    && all(isapprox.([b.model.f_1e for b in cluster.u], ones(2), atol=tol))
  )

  # test case where force is loaded
  cluster, force_on_junction, dt, force_time_dict  = load_from_json(json_path, true)
  force_loaded_ok = (
    cluster isa Cluster
    && cluster.n == 2
    && isapprox(cluster.l, 1.0)
    && cluster.f_model == :force_global
    && all(b -> b isa Bond, cluster.u)
    && all(isapprox.([b.state for b in cluster.u], falses(2), atol=tol))
    && all(isnan.([b.f for b in cluster.u]))
    && all(b -> b.model isa BondModel{Slip}, cluster.u)
    && all(isapprox.([b.model.k_on for b in cluster.u], zeros(2), atol=tol))
    && all(isapprox.([b.model.k_off_0 for b in cluster.u], zeros(2), atol=tol))
    && all(isapprox.([b.model.f_1e for b in cluster.u], ones(2), atol=tol))
    && all(isapprox.(force_on_junction, [1.0, 2.0, 3.0], atol=tol))
    && isapprox(dt, 0.5, atol=tol)
    && length(force_time_dict) == 3
    && haskey(force_time_dict, 0.0)
    && haskey(force_time_dict, 0.5)
    && haskey(force_time_dict, 1.0)
    && isapprox(force_time_dict[0.0][1], 1.0, atol=tol)
    && all(isapprox.(force_time_dict[0.0][2], [0.5, 0.5], atol=tol))
    && isapprox(force_time_dict[0.5][1], 2.0, atol=tol)
    && all(isapprox.(force_time_dict[0.5][2], [1.0, 1.0], atol=tol))
    && isapprox(force_time_dict[1.0][1], 3.0, atol=tol)
    && all(isapprox.(force_time_dict[1.0][2], [1.5, 1.5], atol=tol))
  )

  return no_force_ok && force_loaded_ok

end

@test _check_load_from_json()