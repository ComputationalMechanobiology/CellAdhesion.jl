
println("===============================================")
println("Testing dynamics.jl")
println("===============================================")


function _check_distance()

    l1 = CellAdhesion.distance(BitVector([0,0,1,1,0,1,0,0,1,1,0,1]), 12)
    l2 = CellAdhesion.distance(BitVector([0,0,0,0,0,1,0,0,0,0,0,0]), 12)
    l3 = CellAdhesion.distance(BitVector([0,0,0,0,0,0,0,0,0,0,0,0]), 12)    
    l4 = CellAdhesion.distance(BitVector([1,1,1,0,1,0,1,1,0,0,1,0]), 12)    


    (l1 == [0,0,4,3,0,5,0,0,4,3,0,5]) && (l2 == [0,0,0,0,0,12,0,0,0,0,0,0]) && (l3 == zeros(12)) && (l4 == [3,2,3,0,4,0,3,4,0,0,5,0]) && (typeof(l1) == Vector{CellAdhesionFloat})

end
@test _check_distance()


function _check_force_bonds_global(tol)

    model = Model(Dict("model"=>force_global),Dict(), Dict(), Dict())

    n = 10
    l = 1.0
    v1 = Interface(Bond.([true,true,true,true,true,true,true,true,true,true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)
    v2 = Interface(Bond.([false,true,false,true,false,false,false,false,false,false], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)
    v3 = Interface(Bond.([false,true,true,false,false,false,false,true,false,true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)


    force(v1, model)
    force(v2, model)
    force(v3, model)

    f1 = getfield.(v1.u, :f)
    f2 = getfield.(v2.u, :f)
    f3 = getfield.(v3.u, :f)

    (f1 == repeat([1.0], 10)) && (f2==[0.0, 5.0, 0.0, 5.0, 0.0, 0.0,0.0,0.0,0.0,0.0]) && (f3==[0.0, 2.5, 2.5, 0.0, 0.0, 0.0,0.0, 2.5,0.0, 2.5])

end

@test _check_force_bonds_global(tol)



function _check_force_clusters_global(tol)

    model = Model(Dict("model"=>force_global),Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict())


    int_1 =  interface([2, 3], [1.0, 0.1], 18.0, [true, true], model)
    update_state(int_1)
    force(int_1, model)
    
    f_check_1 = sum(getfield.(int_1.u, :f))
    f_check_2 = 0
    for i = 1:1:int_1.n
        f_check_2 = f_check_2 + sum(getfield.(int_1.u[i].u, :f))
    end

    n = 3
    l = 1
    v1 = Interface(Bond.([false,true,true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
    v2 = Interface(Bond.([false,false,true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
    int_2 = Interface([v1, v2], false, 60.0, false, 2, l)

    update_state(int_2)
    force(int_2, model)
    
    f2_check_1 = sum(getfield.(int_2.u, :f))
    f2_check_2 = 0
    for i = 1:1:int_2.n
        f2_check_2 = f2_check_2 + sum(getfield.(int_2.u[i].u, :f))
    end

    (int_1.f == f_check_1) && (int_1.f == f_check_2) &&  (int_2.f == f2_check_1) && (int_2.f == f2_check_2) 

end

@test _check_force_clusters_global(tol)


function _check_force_bonds_local(tol)

    model = Model(Dict("model"=>force_local),Dict(), Dict(), Dict())

    n = 10
    l = 1.0
    v1 = Interface(Bond.([true,true,true,true,true,true,true,true,true,true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)
    v2 = Interface(Bond.([false,true,false,true,false,false,false,false,false,false], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)
    v3 = Interface(Bond.([false,true,true,false,false,false,false,true,false,true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 10.0, false, n, l)


    force(v1, model)
    force(v2, model)
    force(v3, model)

    f1 = getfield.(v1.u, :f)
    f2 = getfield.(v2.u, :f)
    f3 = getfield.(v3.u, :f)

    (f1 == repeat([1.0], 10)) && (f2==[0.0, 5.0, 0.0, 5.0, 0.0, 0.0,0.0,0.0,0.0,0.0]) && (f3==[0.0, 1.5, 3.0, 0.0, 0.0, 0.0,0.0, 3.5,0.0, 2.0])

end

@test _check_force_bonds_local(tol)




function _check_force_clusters_local(tol)

    model = Model(Dict("model"=>force_global),Dict("model"=>k_on_constant, "k_on_0"=>1.0), Dict("model"=>k_off_slip, "k_off_0"=>0.0, "f_1e"=>1), Dict())



    int_1 =  interface([2, 3], [1.0, 0.1], 18.0, [true, true], model)
    update_state(int_1)
    force(int_1, model)
    
    f_check_1 = sum(getfield.(int_1.u, :f))
    f_check_2 = 0
    for i = 1:1:int_1.n
        f_check_2 = f_check_2 + sum(getfield.(int_1.u[i].u, :f))
    end

    n = 5
    l = 1
    v1 = Interface(Bond.([false,true,true, false, true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
    v2 = Interface(Bond.([false,true,true, true, true], zeros(n), zeros(n), zeros(n), repeat([false], n)), false, 0.0, false, n, l)
    int_2 = Interface([v1, v2], false, 60.0, false, 2, l)

    update_state(int_2)
    force(int_2, model)
    
    f2_check_1 = sum(getfield.(int_2.u, :f))
    f2_check_2 = 0
    for i = 1:1:int_2.n
        f2_check_2 = f2_check_2 + sum(getfield.(int_2.u[i].u, :f))
    end

    (int_1.f == f_check_1) && (int_1.f == f_check_2) &&  (int_2.f == f2_check_1) && (int_2.f == f2_check_2) 

end

@test _check_force_clusters_local(tol)






# function _check_k_off_slip(tol)

#     junction = Interface(false, 4, [0,1,1,1], [0.1, 0.1, 0.1, 0.1], [0.5, 0.5, 0.5, 0.5], [1, 1, 1, 1], BitMatrix(undef,4,0))
#     model = Model(Dict("model"=>k_on_constant, "k_on_0"=>0.2), Dict("model"=>k_off_slip, "k_off_0"=>0.1, "f_1e"=>1), Dict())
#     junction = k_off_slip(junction, model)

#     isapprox(junction.k_off, [0, 0.27182818, 0.27182818, 0.27182818], atol=tol) && (typeof(junction.k_off) == Vector{CellAdhesionFloat})

# end
# @test _check_k_off_slip(tol)


# function _check_k_on_constant(tol)

#     junction = Interface(false, 4, [1,0,0,1], [0.1, 0.1, 0.1, 0.1], [0.5, 0.5, 0.5, 0.5], [1, 1, 1, 1], BitMatrix(undef,4,0))
#     model = Model(Dict("model"=>k_on_constant, "k_on_0"=>0.2), Dict("model"=>k_off_slip, "k_off_0"=>0.1, "f_1e"=>1), Dict())
#     junction = k_on_constant(junction, model)

#     isapprox(junction.k_on, [0,0.2,0.2,0], atol=tol) && (typeof(junction.k_on) == Vector{CellAdhesionFloat})

# end
# @test _check_k_on_constant(tol)

