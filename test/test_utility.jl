
println("===============================================")
println("Testing utility.jl")
println("===============================================")


function _initiate_interface()
    n = 3

    v1 = initiate_interface(n, 1.0, 0.0, true)
    v2 = initiate_interface(n, 1.0, 1.0, true)
    v3 = initiate_interface(n, 1.0, 0.2, true)

    (getproperty.(v1.u[:], :state) == zeros(n)) && (getproperty.(v2.u[:], :state) == ones(n)) && (typeof(getproperty.(v3.u[:], :state)) == BitVector)

end
@test _initiate_interface()



function _initiate_interface()

    v4 = initiate_interface([2, 3], [1.0, 1.1], 0.9, [true, true])

    (length(v4.u)==2) && (length(v4.u[1].u)==3)

end
@test _initiate_interface()






function _update_state()

    int_1 = initiate_interface(5, 1.0, 0.0, true)
    int_2 = initiate_interface(5, 1.0, 1.0, true)
    int_3 = initiate_interface(5, 1.0, 0.5, true)

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)

    int_1.state == false && int_2.state == true && int_3.state == true

end
@test _update_state()


function _update_state()

    int_1 = initiate_interface([2, 3], [1.0, 1.1], 0.0, [true, true])
    int_2 = initiate_interface([2, 3], [1.0, 1.1], 1.0, [true, true])
    int_3 = initiate_interface([2, 3], [1.0, 1.1], 0.5, [true, true])

    update_state(int_1)
    update_state(int_2)
    update_state(int_3)


    (int_1.state == false) && (int_2.state == true) && (int_3.state == true) && (int_1.u[1].state ==false) && (int_1.u[2].state ==false) && (int_2.u[1].state ==true)

end
@test _update_state()


# function _check_distance()

#     l1 = CellAdhesion.distance(BitVector([0,0,1,1,0,1,0,0,1,1,0,1]), 12)
#     l2 = CellAdhesion.distance(BitVector([0,0,0,0,0,1,0,0,0,0,0,0]), 12)
#     l3 = CellAdhesion.distance(BitVector([0,0,0,0,0,0,0,0,0,0,0,0]), 12)    
#     l4 = CellAdhesion.distance(BitVector([1,1,1,0,1,0,1,1,0,0,1,0]), 12)    


#     (l1 == [0,0,4,3,0,5,0,0,4,3,0,5]) && (l2 == [0,0,0,0,0,12,0,0,0,0,0,0]) && (l3 == zeros(12)) && (l4 == [3,2,3,0,4,0,3,4,0,0,5,0]) && (typeof(l1) == Vector{CellAdhesionFloat})

# end
# @test _check_distance()


# function _check_force(tol)

#     model = Model(Dict("model"=>"k_on_constant"), Dict("model"=>"k_off_slip"), Dict("load"=>"global"))
#     f = force(BitVector([0,0,1,1,0,1,0,0,1,1,0,1]), 12, model, convert(CellAdhesionFloat,1))

#     (f == [0,0,2,2,0,2,0,0,2,2,0,2]) && (typeof(f) == Vector{CellAdhesionFloat})

# end

# @test _check_force(tol)




