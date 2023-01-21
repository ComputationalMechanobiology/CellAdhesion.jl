
println("===============================================")
println("Testing utility.jl")
println("===============================================")


function _init_bonds()
    n = 10

    v1 = init_bonds(n, convert(CellAdhesionFloat, 0.0))
    v2 = init_bonds(n, convert(CellAdhesionFloat, 1.0))
    v3 = init_bonds(n, convert(CellAdhesionFloat, 0.2))

    (v1 == zeros(n)) && (v2 == ones(n)) && (typeof(v3) == BitVector)

end
@test _init_bonds()


function _check_state()

    state1 = check_state(BitVector([0,0,0,0,0]))
    state2 = check_state(BitVector([0,0,1,0,0]))

    state1 == true && state2 == false

end
@test _check_state()



function _check_distance()

    l1 = CellAdhesion.distance(BitVector([0,0,1,1,0,1,0,0,1,1,0,1]), 12)
    l2 = CellAdhesion.distance(BitVector([0,0,0,0,0,1,0,0,0,0,0,0]), 12)
    l3 = CellAdhesion.distance(BitVector([0,0,0,0,0,0,0,0,0,0,0,0]), 12)    
    l4 = CellAdhesion.distance(BitVector([1,1,1,0,1,0,1,1,0,0,1,0]), 12)    


    (l1 == [0,0,4,3,0,5,0,0,4,3,0,5]) && (l2 == [0,0,0,0,0,12,0,0,0,0,0,0]) && (l3 == zeros(12)) && (l4 == [3,2,3,0,4,0,3,4,0,0,5,0]) && (typeof(l1) == Vector{CellAdhesionFloat})

end
@test _check_distance()


function _check_force(tol)

    model = Model(Dict("model"=>"k_on_constant"), Dict("model"=>"k_off_slip"), Dict("load"=>"global"))
    f = force(BitVector([0,0,1,1,0,1,0,0,1,1,0,1]), 12, model, convert(CellAdhesionFloat,1))

    (f == [0,0,2,2,0,2,0,0,2,2,0,2]) && (typeof(f) == Vector{CellAdhesionFloat})

end

@test _check_force(tol)




