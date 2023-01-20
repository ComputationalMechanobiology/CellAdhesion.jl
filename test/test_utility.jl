
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

    junction1 = Interface(false, 5, [0,0,0,0,0], [], [], [], BitMatrix(undef,5,0))
    junction1 = check_state(junction1)

    junction2 = Interface(true, 5, [0,0,1,0,0], [], [], [], BitMatrix(undef,5,0))
    junction2 = check_state(junction2)

    junction1.state == true && junction2.state == false

end
@test _check_state()



function _check_distance()

    junction1 = Interface(false, 12, [0,0,1,1,0,1,0,0,1,1,0,1], [], [], [], BitMatrix(undef,12,0))
    l1 = CellAdhesion.distance(junction1)

    junction2 = Interface(false, 12, [0,0,0,0,0,1,0,0,0,0,0,0], [], [], [], BitMatrix(undef,12,0))
    l2 = CellAdhesion.distance(junction2)

    junction3 = Interface(false, 12, [0,0,0,0,0,0,0,0,0,0,0,0], [], [], [], BitMatrix(undef,12,0))
    l3 = CellAdhesion.distance(junction3)    

    junction4 = Interface(false, 12, [1,1,1,0,1,0,1,1,0,0,1,0], [], [], [], BitMatrix(undef,12,0))
    l4 = CellAdhesion.distance(junction4)    


    (l1 == [0,0,4,3,0,5,0,0,4,3,0,5]) && (l2 == [0,0,0,0,0,12,0,0,0,0,0,0]) && (l3 == zeros(12)) && (l4 == [3,2,3,0,4,0,3,4,0,0,5,0]) && (typeof(l1) == Vector{CellAdhesionFloat})

end
@test _check_distance()


function _check_force(tol)

    junction = Interface(false, 12, [0,0,1,1,0,1,0,0,1,1,0,1], [], [], [], BitMatrix(undef,12,0))
    param = Dict("load"=>"global")
    junction = force(junction, param, convert(CellAdhesionFloat,1))

    (junction.f == [0,0,2,2,0,2,0,0,2,2,0,2]) && (typeof(junction.f) == Vector{CellAdhesionFloat})

end

@test _check_force(tol)




