
export initiate_interface, update_state, force, force_global, force_local



function initiate_interface(n::Union{Int, Vector{Int}}, l::Union{Float64, Vector{Float64}}, K::Float64, F::Float64, history::Union{Bool, Vector{Bool}})

  check_length = length(n)
  @assert length(l) == check_length "Length missmatch initialisation"
  @assert length(K) == 1 "Length missmatch initialisation"
  @assert length(history) == check_length "Length missmatch initialisation"

  # Make variable types consistent with CellAdhesionFloat
  if check_length == 1
    n = convert(Int32, n)
    l = convert(CellAdhesionFloat, l)
  else
    n = convert(Vector{Int32}, n)
    l = convert(Vector{CellAdhesionFloat}, l)
  end
    K = convert(CellAdhesionFloat, K)
    F = convert(CellAdhesionFloat, F)

  # If the junction has only one cluster
  if check_length ==1

    x = Interface(init_bonds(n, K, history), false, F, false, n, l)
  
  else  # If the junction has multiple clusters

    x = Interface(Vector{Interface}(undef,0), false, F, false, n[1], l[1])
    for i = 1:1:n[1]
      push!(x.u, Interface(init_bonds(n[2], K, history[2]), false, 0.0, history[1], n[2], l[2]))
    end

  end
  
  return x

end

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
distance(v::BitVector, n::Integer)

  Compute the distance of each from its two closest closed link on each side. 
  Used to compute a "local" load distributions across closed bonds.
  Periodic bondary conditions for the edges.

  Input parameters:
    - v: vector with the state of each single bond
    - n: number of bonds in the junction
  Output parameters:
    - l: vector of CellAdhesionFloat containing the distance for each closed bond
"""

function distance(v::BitVector, n::Integer)

  @assert !isempty(v) "Bond state vector in Interface is empty"

  l = zeros(n);
  if sum(v) == 0
       l .= zeros(n);
  elseif sum(v) ==1
      temp = findall(x->x==1, v)
      l[temp[1]] = n;        
  elseif sum(v) == 2
      temp = findall(x->x==1, v)
      l[temp[1]] = n;
      l[temp[2]] = n;    
  else
      temp = findall(x->x==1, v)
      gaps = diff(temp);
      dist = zeros(size(temp));
      dist[2:end-1] = gaps[2:end] .+ gaps[1:end-1];
      dist[1] = n - temp[end] + temp[2];
      dist[end] = n - temp[end-1] + temp[1];
      l[temp] .= dist;
  end

  return convert(Vector{CellAdhesionFloat}, l)

end


"""
force(v::Interface, model::Model, s::CellAdhesionFloat)

  Compute the force on each link. 
  Two options are available:
   - global: equal distribution across all closed bonds
   - local: inhomogeneous distribution of load accounting fro the distance to its nearest closed neighbor on both sides. 

   This is defined in the model.f variable

  Input parameters:
    - v: Interface structure
    - model: Model structure containing the type of load distribution: "local" or "global" (Model.f)
  Output parameters:
    - Updated Interface with force applied to each link
"""


function force(v::Interface, force_type)

  @assert v.f>=0 "Applied stress to junction must be positive or equal to zero"

  force_type(v, v.f)

  if typeof(v.u)== Vector{Interface}
    
    for i = 1:1:v.n
      force_type(v.u[i], v.u[i].f)
    end
  
  end

end



function force_global(v::Interface, f::CellAdhesionFloat)

  interface_v = getfield.(v.u, :state);
  update_f = interface_v .* f./sum(interface_v)
  setproperty!(v, :f, convert(Vector{CellAdhesionFloat},update_f))

end


function force_local(v::Interface, f::CellAdhesionFloat)

  interface_v = getfield.(v.u, :state);
  l = distance(interface_v, v.n)
  update_f = l ./ sum(l) .*f;
  setproperty!(v, :f, convert(Vector{CellAdhesionFloat},update_f))

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


