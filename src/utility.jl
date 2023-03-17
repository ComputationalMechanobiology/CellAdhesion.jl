
export update_state, KineticMonteCarlo_unit, KineticMonteCarlo




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
Base.setproperty!    NOT TESTED YET!

Definition to update Interface struct fields
"""


function Base.setproperty!(x::Interface, s::Symbol, new_x::Vector{CellAdhesionFloat})

  for i = 1:1:x.n
    setfield!(x.u[i], s, new_x[i])
  end

end




function KineticMonteCarlo(v::Interface, model::Model)

  #@assert (model.param["dt"]>0) "dt not valid"
  
  if typeof(v.u)==Vector{Bond}
    KineticMonteCarlo_unit(v, model.param["dt"])
  else

    for i=1:1:v.n
      KineticMonteCarlo_unit(v.u[i], model.param["dt"])
    end

  end
    
end




"""
KineticMonteCarlo(v::Interface, dt::CellAdhesionFloat)

"""

function KineticMonteCarlo_unit(v::Interface, dt::CellAdhesionFloat)

  random = rand(v.n);
  v_temp = repeat([false],v.n);

  bond_events = findall(x->x==1, (getfield.(v.u, :k_on) .* dt .>random));
  v_temp[bond_events] .= true;
  unbond_events = findall(x->x==1, (getfield.(v.u, :k_off).*dt .>random));
  v_temp[unbond_events] .= false;

  for i = 1:1:v.n                               # if I use setproperty! it doesn't work! Ask for help!
    setfield!(v.u[i], :state, v_temp[i])
  end

end


