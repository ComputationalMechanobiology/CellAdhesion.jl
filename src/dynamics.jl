
export force, k_on_constant, k_off_slip, k_rate_junction




#------------------ FORCE -------------------------------------

"""
force(v::Interface, model::Model)

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


function force(v::Cluster)

  @assert v.f>=0 "Applied stress to Cluster must be positive or equal to zero"

  getfield(CellAdhesion, v.f_model)(v)      # If we define v.f_model as string then => Symbol(v.f_model)

end


function force(v::Interface)

  @assert v.f>=0 "Applied stress to Interface must be positive or equal to zero"

  getfield(CellAdhesion, Symbol(v.f_model))(v)
  
  for i = 1:1:v.n
    getfield(CellAdhesion, Symbol(v.u[i].f_model))(v.u[i])
  end
  
end

"""
force_global
Computer force distribution by equally dividing the force within the closed bonds
"""

function force_global(v::Union{Cluster, Interface})

  interface_v = getfield.(v.u, :state);
  update_f = interface_v .* v.f./sum(interface_v)
  setproperty!(v, :f, convert(Vector{CellAdhesionFloat},update_f))

end

"""
force_local
Computer force distribution by accounting for the distance of each link from its two closest neighbours
"""


function force_local(v::Union{Cluster, Interface})

  interface_v = getfield.(v.u, :state);
  l = distance(interface_v, v.n)
  update_f = l ./ sum(l) .*v.f;
  setproperty!(v, :f, convert(Vector{CellAdhesionFloat},update_f))

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



#------------------ K RATE -------------------------------------

function k_rate_junction(v::Cluster)

  k_on = getindex(getfield(v.u[1], :model))
  interim = getfield(k_on, :k_on)
  model_k_on = interim[:model]

  k_off = getindex(getfield(v.u[1], :model))
  interim = getfield(k_off, :k_off)
  model_k_off = interim[:model]
  
  model_k_on(v)
  model_k_off(v)
 
end


function k_rate_junction(v::Interface)

  k_on = getindex(getfield(v.u[1].u[1], :model))
  interim = getfield(k_on, :k_on)
  model_k_on = interim[:model]

  k_off = getindex(getfield(v.u[1].u[1], :model))
  interim = getfield(k_off, :k_off)
  model_k_off = interim[:model]
  
  for i=1:1:v.n
    model_k_on(v.u[i])
    model_k_off(v.u[i])
  end

end


"""
k_off_slip(v::Union{Cluster, Interface})

  Compute the probability of unbinding of each closed bonds within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - v: Interface variable
    - model: Model varible containing the unbinding parameters (k_off_0 = rate of unbinding with no force applied, 
                                                                f_1e = threshold force exponential decay)
  Output parameters:
    - k_off: Unbinding probability for each bond
"""

function k_off_slip(v::Union{Cluster, Interface})

  k_off = getindex.(getfield.(v.u, :model))
  interim = getfield.(k_off, :k_off)
  model_vect, k_off_0_vect, f_1e_vect = [getproperty.(interim, i) for i in (:model, :k_off_0, :f_1e)]
  update_k_off = k_off_0_vect .* exp.(getfield.(v.u, :f) ./ f_1e_vect) .* getfield.(v.u, :state);
  setproperty!(v, :k_off, convert(Vector{CellAdhesionFloat},update_k_off))

end


"""
k_on_constant(v::Union{Cluster, Interface})

  Compute the probability of binding of each open bond within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - v: Interface variable
    - model: Model varible containing the binding parameters (k_on_0 = rate of binding with no force applied)
  Output parameters:
    - junction: binding probability for each bond
"""

function k_on_constant(v::Union{Cluster, Interface})

  k_on = getindex.(getfield.(v.u, :model))
  interim = getfield.(k_on, :k_on)
  # # keys_k_on = keys(interim[1])
  # # print(keys_k_on)
  model_vect, k_on_0_vect = [getproperty.(interim, i) for i in (:model, :k_on_0)]     #X, Y, Z = [getindex.(interim, i) for i in 1:2]  #number of keys
  update_k_on = k_on_0_vect .* (ones(v.n) .- getfield.(v.u, :state));
  setproperty!(v, :k_on, convert(Vector{CellAdhesionFloat},update_k_on))

end


# Add here new k_on or k_off functions models!

list_k_models = [k_off_slip, 
                 k_on_constant]


