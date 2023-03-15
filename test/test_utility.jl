
println("===============================================")
println("Testing utility.jl")
println("===============================================")


function _initiate_interface()
    n = 3

    v1 = initiate_interface(n, 1.0, 0.0, 15.0, true)
    v2 = initiate_interface(n, 1.0, 1.0, 15.0, true)
    v3 = initiate_interface(n, 1.0, 0.2, 15.0, true)

    (getproperty.(v1.u[:], :state) == zeros(n)) && (getproperty.(v2.u[:], :state) == ones(n)) && (typeof(getproperty.(v3.u[:], :state)) == BitVector) && (v1.f == 15.0)

end
@test _initiate_interface()



function _initiate_interface()

    v4 = initiate_interface([2, 3], [1.0, 1.1], 0.9, 15.0, [true, true])

    (length(v4.u)==2) && (length(v4.u[1].u)==3) && (v4.f == 15.0)

end
@test _initiate_interface()






function _update_state_bonds()

    int_1 = initiate_interface(5, 1.0, 0.0, 15.0, true)
    int_2 = initiate_interface(5, 1.0, 1.0, 15.0, true)
    int_3 = initiate_interface(5, 1.0, 0.5, 15.0, true)

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)

    int_1.state == false && int_2.state == true && int_3.state == true

end
@test _update_state_bonds()


function _update_state_clusters()

    int_1 = initiate_interface([2, 3], [1.0, 1.1], 0.0, 15.0, [true, true])
    int_2 = initiate_interface([2, 3], [1.0, 1.1], 1.0, 15.0,[true, true])
    int_3 = initiate_interface([2, 3], [1.0, 1.1], 0.5, 15.0, [true, true])

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)


    (int_1.state == false) && (int_2.state == true) && (int_3.state == true) && (int_1.u[1].state ==false) && (int_1.u[2].state ==false) && (int_2.u[1].state ==true)

end
@test _update_state_clusters()


function _check_distance()

    l1 = CellAdhesion.distance(BitVector([0,0,1,1,0,1,0,0,1,1,0,1]), 12)
    l2 = CellAdhesion.distance(BitVector([0,0,0,0,0,1,0,0,0,0,0,0]), 12)
    l3 = CellAdhesion.distance(BitVector([0,0,0,0,0,0,0,0,0,0,0,0]), 12)    
    l4 = CellAdhesion.distance(BitVector([1,1,1,0,1,0,1,1,0,0,1,0]), 12)    


    (l1 == [0,0,4,3,0,5,0,0,4,3,0,5]) && (l2 == [0,0,0,0,0,12,0,0,0,0,0,0]) && (l3 == zeros(12)) && (l4 == [3,2,3,0,4,0,3,4,0,0,5,0]) && (typeof(l1) == Vector{CellAdhesionFloat})

end
@test _check_distance()


function _check_force_bonds_global(tol)

    model = Model(Dict("model"=>"k_on_constant"), Dict("model"=>"k_off_slip"), Dict("load"=>"global"))

    n = 10
    l = 1.0
    v1 = initiate_interface(10, 1.0, 1.0, 10.0, true)
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



function _check_force_clusters(tol)

    model = Model(Dict("model"=>"k_on_constant"), Dict("model"=>"k_off_slip"), Dict("load"=>"global"))


    int_1 =  initiate_interface([2, 3], [1.0, 1.1], 1.0, 18.0, [true, true])
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




