
println("===============================================")
println("Testing dynamics.jl")
println("===============================================")


function _check_k_on(tol)

  model = SlipBondModel((k_on_0 = 0.2,), (k_off_0 = 0.8, f_1e = 1.0))
  k_on_bond = CellAdhesion.k_on(model)

  k_on_bond == 0.2

end

_check_k_on(tol)


function _check_k_off(tol)

  model = SlipBondModel((k_on_0 = 0.2,), (k_off_0 = 0.8, f_1e = 1.0))
  k_off_bond = CellAdhesion.k_off(model, convert(CellAdhesionFloat, 2.0))

  k_off_bond == 0.8*exp(2)

end

_check_k_off(tol)



function _check_distance()

    l1 = CellAdhesion.distance(BitVector([false,false,true, true, false,true, false,false,true, true, false, true]), 12)
    l2 = CellAdhesion.distance(BitVector([false,false,false,false,false,true, false,false,false,false,false,false]), 12)
    l3 = CellAdhesion.distance(BitVector([false,false,false,false,false,false,false,false,false,false,false,false]), 12)    
    l4 = CellAdhesion.distance(BitVector([true, true, true, false,true, false,true, true, false,false,true, false]), 12)    


    ((l1 == [0,0,4,3,0,5,0,0,4,3,0,5]) 
      && (l2 == [0,0,0,0,0,12,0,0,0,0,0,0]) 
      && (l3 == zeros(12)) 
      && (l4 == [3,2,3,0,4,0,3,4,0,0,5,0]) 
      && (typeof(l1) == Vector{CellAdhesionFloat}))

end
@test _check_distance()


function _check_force_bonds_global(tol)

    force_string = :force_global
  
    model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  
    n = convert(CellAdhesionInt, 10)
    l = convert(CellAdhesionFloat, 1.0)
    v1 = Cluster(Bond.([true,true,true,true,true,true,true,true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
    v2 = Cluster(Bond.([false,true,false,true,false,false,false,false,false,false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
    v3 = Cluster(Bond.([false,true,true,false,false,false,false,true,false,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)

    F = convert(CellAdhesionFloat, 10.0)
    setforce!(v1, F)
    setforce!(v2, F)
    setforce!(v3, F)

    f1 = getfield.(v1.u, :f)
    f2 = getfield.(v2.u, :f)
    f3 = getfield.(v3.u, :f)

    ((f1 == repeat([1.0], 10)) 
      && (f2==[0.0, 5.0, 0.0, 5.0, 0.0, 0.0,0.0,0.0,0.0,0.0]) 
      && (f3==[0.0, 2.5, 2.5, 0.0, 0.0, 0.0,0.0, 2.5,0.0, 2.5]))

end

@test _check_force_bonds_global(tol)



function _check_force_clusters_global(tol)

  model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  n = convert(CellAdhesionInt, 3)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 60.0)

  force_string = :force_global
  v1 = Cluster(Bond.([true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)

  setforce!(int_1, F)

  
  f_check_1 = sum(getfield.(int_1.u, :f))
  f_check_2 = 0
  for i = 1:1:int_1.n
      f_check_2 = f_check_2 + sum(getfield.(int_1.u[i].u, :f))
  end

  v1 = Cluster(Bond.([false,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([false,false,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_2 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  
  setforce!(int_2, F)
  
  f2_check_1 = sum(getfield.(int_2.u, :f))
  f2_check_2 = 0
  for i = 1:1:int_2.n
      f2_check_2 = f2_check_2 + sum(getfield.(int_2.u[i].u, :f))
  end


  ((int_1.f == f_check_1) 
    && (int_1.f == f_check_2) 
    && (int_2.f == f2_check_1) 
    && (int_2.f == f2_check_2)) 

end

@test _check_force_clusters_global(tol)


function _check_force_clusters2levels_global(tol)

  model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  n = convert(CellAdhesionInt, 4)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 60.0)

  force_string = :force_global
  v1 = Cluster(Bond.([true,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([false,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v3 = Cluster(Bond.([true,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v4 = Cluster(Bond.([true,true,true, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v5 = Cluster(Bond.([false,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v6 = Cluster(Bond.([true,true,true, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)

  c1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  c2 = Cluster([v3, v4], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  c3 = Cluster([v5, v6], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  
  int_1 = Cluster([c1, c2, c3], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 3), l)

  setforce!(int_1, F)

  f_check_1 = sum(getfield.(int_1.u, :f))
  f_check_2 = 0
  for i = 1:1:int_1.n
      f_check_2 = f_check_2 + sum(getfield.(int_1.u[i].u, :f))
  end


  ((int_1.f == f_check_1) 
    && (int_1.f == f_check_2)) 

end

@test _check_force_clusters2levels_global(tol)



function _check_force_bonds_local(tol)


  force_string = :force_local
  
  model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))

  n = convert(CellAdhesionInt, 10)
  l = convert(CellAdhesionFloat, 1.0)
  v1 = Cluster(Bond.([true,true,true,true,true,true,true,true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([false,true,false,true,false,false,false,false,false,false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v3 = Cluster(Bond.([false,true,true,false,false,false,false,true,false,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)

  F = convert(CellAdhesionFloat, 10.0)
  setforce!(v1, F)
  setforce!(v2, F)
  setforce!(v3, F)

  f1 = getfield.(v1.u, :f)
  f2 = getfield.(v2.u, :f)
  f3 = getfield.(v3.u, :f)

  ((f1 == repeat([1.0], 10)) 
    && (f2==[0.0, 5.0, 0.0, 5.0, 0.0, 0.0,0.0,0.0,0.0,0.0]) 
    && (f3==[0.0, 1.5, 3.0, 0.0, 0.0, 0.0,0.0, 3.5,0.0, 2.0]))

end

@test _check_force_bonds_local(tol)


function _check_force_clusters_local(tol)

  model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  n = convert(CellAdhesionInt, 3)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 18.0)

  force_string = :force_local
  v1 = Cluster(Bond.([true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)

  setforce!(int_1, F)

  
  f_check_1 = sum(getfield.(int_1.u, :f))
  f_check_2 = 0
  for i = 1:1:int_1.n
      f_check_2 = f_check_2 + sum(getfield.(int_1.u[i].u, :f))
  end

  F = convert(CellAdhesionFloat, 60.0)
  v1 = Cluster(Bond.([false,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([false,false,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  int_2 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  
  setforce!(int_2, F)
  
  f2_check_1 = sum(getfield.(int_2.u, :f))
  f2_check_2 = 0
  for i = 1:1:int_2.n
      f2_check_2 = f2_check_2 + sum(getfield.(int_2.u[i].u, :f))
  end


  ((int_1.f == f_check_1) 
    && (int_1.f == f_check_2) 
    && (int_2.f == f2_check_1) 
    && (int_2.f == f2_check_2)) 

end

@test _check_force_clusters_local(tol)



function _check_force_clusters2levels_local(tol)

  model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
  n = convert(CellAdhesionInt, 4)
  l = convert(CellAdhesionFloat, 1.0)
  F = convert(CellAdhesionFloat, 60.0)

  force_string = :force_local
  v1 = Cluster(Bond.([true,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v2 = Cluster(Bond.([false,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v3 = Cluster(Bond.([true,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v4 = Cluster(Bond.([true,true,true, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v5 = Cluster(Bond.([false,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
  v6 = Cluster(Bond.([true,true,true, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)

  c1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  c2 = Cluster([v3, v4], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  c3 = Cluster([v5, v6], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  
  int_1 = Cluster([c1, c2, c3], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 3), l)

  setforce!(int_1, F)

  
  f_check_1 = sum(getfield.(int_1.u, :f))
  f_check_2 = 0
  for i = 1:1:int_1.n
      f_check_2 = f_check_2 + sum(getfield.(int_1.u[i].u, :f))
  end


  ((int_1.f == f_check_1) 
    && (int_1.f == f_check_2)) 

end

@test _check_force_clusters2levels_local(tol)


function _check_force_clusters2levels_local_test2(tol)

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

  
  f_check_1 = sum(getfield.(int_1.u, :f))
  f_check_2 = 0
  for i = 1:1:int_1.n
      f_check_2 = f_check_2 + sum(getfield.(int_1.u[i].u, :f))
  end


  ((int_1.f == f_check_1) 
    && (int_1.f == f_check_2)) 

end

@test _check_force_clusters2levels_local_test2(tol)





