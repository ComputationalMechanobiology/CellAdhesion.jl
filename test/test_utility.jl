
println("===============================================")
println("Testing utility.jl")
println("===============================================")


function _check_SlipBondModel()

    k_on_params = (k_on_0 = 0.2,)
    k_off_params = (k_off_0 = 0.8, f_1e = 1.0)

    model = SlipBondModel(k_on_params, k_off_params)

    typeof(model) == SlipBondModel

end

@test _check_SlipBondModel()



function _check_print()
    model1 = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
    model2 = SlipBondModel((k_on_0=0.2,), (k_off_0=0.8, f_1e=1))
    n = convert(CellAdhesionInt, 3)
    l = convert(CellAdhesionFloat, 1.0)
    F = convert(CellAdhesionFloat, 60.0)
  
    force_string = :force_global
    v1 = Cluster(Bond.([true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model1], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
    v2 = Cluster(Bond.([true,true,true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model2], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
    int_1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
  

    #print_cluster(int_1)

    1==1

end

@test _check_print()



function _check_state()

    model = SlipBondModel((k_on_0=1.0,), (k_off_0=0.0, f_1e=1))
    n = convert(CellAdhesionInt, 4)
    l = convert(CellAdhesionFloat, 1.0)
    F = convert(CellAdhesionFloat, 60.0)
  
    force_string = :force_global
    v1 = Cluster(Bond.([true,false,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
    v2 = Cluster(Bond.([false,false,false, false], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), true, convert(CellAdhesionFloat, 0.0), force_string, n, l)
    v3 = Cluster(Bond.([true,true,true, true], convert(Vector{CellAdhesionFloat}, zeros(n)), repeat([model], n)), false, convert(CellAdhesionFloat, 0.0), force_string, n, l)
    c1 = Cluster([v1, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
    c2 = Cluster([v1, v3], false, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
    c3 = Cluster([v2, v2], true, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 2), l)
    
    int_1 = Cluster([c1, c2, c3], false, convert(CellAdhesionFloat, 0.0), force_string, convert(CellAdhesionInt, 3), l)

    CellAdhesion.state!(int_1)

    ((int_1.state == true) 
     && (int_1.u[1].state == true)
     && (int_1.u[2].state == true)
     && (int_1.u[3].state == false)
     && (int_1.u[1].u[1].state == true)
     && (int_1.u[1].u[2].state == false)
     && (int_1.u[2].u[1].state == true)
     && (int_1.u[2].u[2].state == true)
     && (int_1.u[3].u[1].state == false)
     && (int_1.u[3].u[2].state == false))

end

@test _check_state()