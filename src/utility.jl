
export update_state




"""
init_bonds(n::Integer, K::CellAdhesionFloat)

  Generate the initial state of each bond within the cell junction.

  Input parameters:
    - n: number of bonds in the junction
    - K: probability to start in a closed state
  Output parameters:
    - v: state vecotr of individual bonds (0 = open, 1 = closed)
"""

function init_bonds(n::Integer,K::CellAdhesionFloat, history::Bool)

    @assert n>0 "Number of bonds cannot be negative or equal to 0"
    @assert (K>=0) && (K<=1) "Initialisation probability must be 0<=K<=1"

    v = isless.(rand(n),K)

    return Bond.(v, zeros(n), zeros(n), zeros(n), repeat([history], n))

    return v

end


"""
update_unit_state(v::Interface)

  update unit junction state - function not exported, called by update_state()

"""

function update_unit_state(v::Interface)

  # Get the state value for each bond
  interface_v = getfield.(v.u, :state);

  # If the sum of the state values is 0, the junction is broken 
  sum_v = sum(interface_v);
  state = isequal(sum_v,0);

  # Update the state value of the junction
  setfield!(v, :state, !state)

  return !state
  

end


"""
update_state(v::Interface)

  update if a junction is broken or still viable.

  Input parameters:
    - v: Interface object 
  Output parameters:
    - v with updated state field
"""



function update_state(v::Interface)

  @assert !isempty(v.u) "Unit vector in Interface is empty"

  if typeof(v.u)== Vector{Interface}

    for i = 1:1:v.n
      update_unit_state(v.u[i])
    end

  end
    
  update_unit_state(v)


end

"""
Base.setproperty! 

Definition to update Interface struct fields
"""


function Base.setproperty!(x::Interface, s::Symbol, new_x::Vector{CellAdhesionFloat})

  for i = 1:1:x.n
    setfield!(x.u[i], s, new_x[i])
  end

end






"""
KineticMonteCarlo(v::BitVector, n::Integer, k_on::Vector{CellAdhesionFloat}, k_off::Vector{CellAdhesionFloat}, model::Model)-NOT TESTED YET!!!!!

"""

function KineticMonteCarlo(v::BitVector, n::Integer, k_on::Vector{CellAdhesionFloat}, k_off::Vector{CellAdhesionFloat}, model::Model)

  random = rand(n);
  v_temp = copy(v);

  bond_events = findall(x->x==1, ((k_on.*model.param["dt"]) .>random));
  v_temp[bond_events] .= 1;
  unbond_events = findall(x->x==1, ((k_off.*model.param["dt"]) .>random));
  v_temp[unbond_events] .= 0;


  return v_temp


end


