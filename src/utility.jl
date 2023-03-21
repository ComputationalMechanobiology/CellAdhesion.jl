
export model_init, update_state, KineticMonteCarlo_unit, KineticMonteCarlo



function model_init(force_dict::Dict, k_on_dict::Dict, k_off_dict::Dict, param_dict::Dict)

  return Ref(Model(force_dict, k_on_dict, k_off_dict, param_dict))

end


"""
init_bonds(n::Integer, l::CellAdhesionFloat, F::CellAdhesionFloat, model::Base.RefValue{Model})

  Generate the initial state of each bond within the cell junction.

  Input parameters:
    - n: number of bonds in the junction
    - 
  Output parameters:
    - Interface struct
"""

function init_bonds(n::CellAdhesionInt, l::CellAdhesionFloat, F::CellAdhesionFloat, model::Base.RefValue{Model})

  K = model[].k_on["k_on_0"] / (model[].k_on["k_on_0"] + model[].k_off["k_off_0"])
  v = isless.(rand(n),K)

  return Interface(Bond.(v, zeros(n), zeros(n), zeros(n), repeat([model],n)), false, F, n, l)

end


function init_bonds(n::Vector{CellAdhesionInt}, l::Vector{CellAdhesionFloat}, F::CellAdhesionFloat, model::Base.RefValue{Model})

  v = Interface(Vector{Interface}(undef,n[1]), false, F, n[1], l[1])

  for i = 1:1:n[1]
    v.u[i] = init_bonds(n[2], l[2], F, model)
  end

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
  v_temp = getfield.(v.u, :state);

  bond_events = findall(x->x==1, (getfield.(v.u, :k_on) .* dt .>random));
  v_temp[bond_events] .= true;
  unbond_events = findall(x->x==1, (getfield.(v.u, :k_off).*dt .>random));
  v_temp[unbond_events] .= false;

  for i = 1:1:v.n                               # if I use setproperty! it doesn't work! Ask for help!
    setfield!(v.u[i], :state, v_temp[i])
  end

end


