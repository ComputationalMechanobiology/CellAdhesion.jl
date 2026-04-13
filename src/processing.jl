export update!, runcluster, Cluster, Bond, bond_state_force, save_cluster_to_json, load_from_json

"""
update!(v, dt)

Updates the state of an interface (that can be made of bonds or it can be a hierarchical structure).

Input parameters:
 - v: interface (this can be a cluster or a bond)
 - dt: time step of the simulation

"""
function update!(v::Bond, dt::CellAdhesionFloat)

  random = rand();
  v.state ? k = k_off(v.model, v.f) : k = k_on(v.model)
  (k*dt>random) ? (v.state = !(v.state)) : nothing

end


function update!(v::Cluster, dt::CellAdhesionFloat)

    for i = 1:1:v.n
      k = v.u[i]
      update!(k, dt)
    end

    #Get the state value of the cluster
    clusterstate = false
    for i = 1:1:v.n
      clusterstate |= v.u[i].state
    end

    # Update the state value of the junction
    setfield!(v, :state, clusterstate)

end



"""
runcluster(v, force, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false)

Simulates a junction subjected to an external force using a Montecarlo algorithm. 

Input paramters:
  - v: structure of type Cluster
  - force: it can either be an constant number (if the junction is subjected to a constant force), or a vector (if the junction is subjected to a varing force)
  - dt: time step for the simulation
  - max_steps: maximum number of iterations if the junction doesn't break
  
Output paramters:
  - state of the whole Cluster 
  - force at which it breaks
  - time at which it breaks
  - number of steps after which it breaks

"""
function runcluster(v::Cluster, force::Float64, dt::Float64; state_check::Bool = true, max_steps::Integer = 1000, verbose::Bool = false)

  step = 0
  force = convert(CellAdhesionFloat,force)
  dt = convert(CellAdhesionFloat, dt)

  if state_check == true
    while (step <= max_steps) && (v.state == true)
      step = step + 1
      setforce!(v, force)
      update!(v, dt)
    end
  else
    while (step <= max_steps)
      step = step + 1
      setforce!(v, force)
      update!(v, dt)
    end
  end


  if verbose == true
      if v.state == false
          print("Junction broken")
      elseif step > max_steps
          print("Maximum number of iterations reached")
      end
  end

  return v.state, force, dt*step, step

end

function runcluster(v::Cluster, force::Vector{Float64}, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false)

  # Arbitrary force history applied to the junction
  n = length(force)

  if max_steps > n
         @warn max_steps<=n "Maximum number of steps exceed force vector length"
  	 max_steps = n
         print("\n Maximum number of steps = ", max_steps, "\n")
  end


  step = 1

  force = convert(Vector{CellAdhesionFloat},force)
  dt = convert(CellAdhesionFloat, dt)

  while (step <= max_steps) && (v.state == true)
      F = force[step]
      setforce!(v, F)
      update!(v, dt)
      step = step + 1
  end

  if verbose == true
      if v.state == false
          print("Junction broken")
      elseif step > max_steps
          print("Maximum number of iterations reached")
      end
  end

  return v.state, force[step-1], dt*(step-1), (step-1)


end


function runcluster(v::Cluster, force::Vector{Float64}, dt::Float64, json_file_name::String; max_steps::Integer = 1000, verbose::Bool = false)
  # Arbitrary force history applied to the junction
  n = length(force)
  if max_steps > n
         @warn max_steps<=n "Maximum number of steps exceed force vector length"
  	 max_steps = n
         print("\n Maximum number of steps = ", max_steps, "\n")
  end
  
  step = 1
  states = Vector{AbstractDict}()
  force = convert(Vector{CellAdhesionFloat},force)
  dt = convert(CellAdhesionFloat, dt)
  while (step <= max_steps) && (v.state == true)
      F = force[step]
      setforce!(v, F)
      push!(states, bond_state_force(v; output = :nested, time = (step-1)*dt))
      update!(v, dt)
      step = step + 1
  end
  if v.state == false
    # add the state after the final update!
    F = force[step-1]
    setforce!(v, F)
    push!(states, bond_state_force(v; output = :nested, time = (step-1)*dt))
  end
  if verbose == true
    if v.state == false
        print("Junction broken")
    elseif step > max_steps
        print("Maximum number of iterations reached")
    end
  end

  # save results to json file
  combined_dict = Dict()
  # save each time point's force information
  for state in states
      for (time, force_tree) in state
          combined_dict[time] = force_tree
      end
  end
  cluster_dict = _save_cluster_to_dict(v)   # save cluster structure information
  combined_dict["cluster"] = cluster_dict
  combined_dict["force"] = force
  combined_dict["dt"] = dt

  # save to json file
  open(json_file_name, "w") do io
    JSON.json(io, combined_dict; allownan=true, nan="nan")
  end

  return v.state, force[step-1], dt*(step-1), (step-1)
end




"""
Cluster(n, l, model, f_model)

Initiates a Cluster structure

Input paramters:
 - n: if it is a single number it represents the number of bonds in the junction, if it is a vector it creates a hierarchical structure 
      where each element in the vector represents the number of subunits per hierarchical level. 
 - l: if it is a single number it represents the distance between bonds, if it is a vector each element is the distance between the subunits
 - model: type of BondModel to compute the probability of binding and unbinding
 - f_model: "global" or "local". It defines the method used to redistribute the force across subunits (Bonds or other Clusters if it is a hierarchical structure)

Output paramters:
 - Cluster data structure
 """
 
function Cluster(n::CellAdhesionInt, l::CellAdhesionFloat, model::T, f_model::Symbol) where T<:BondModel

  u = Vector{Bond{T}}(undef, n)

  for i = 1:1:n
    u[i] = Bond(model)
  end
  x = Cluster(u, false, convert(CellAdhesionFloat,NaN), f_model, n, l)
  
  state!(x)

  return x

end

function Cluster(n::Vector{CellAdhesionInt}, l::Vector{CellAdhesionFloat}, model::T, f_model::Vector{Symbol}) where T<:BondModel

  u = Vector{Cluster}(undef, n[1])

  for i = 1:1:n[1]
    if length(n) == 2
      u[i] = Cluster(n[2],l[2], model, f_model[2])
    else
      u[i] = Cluster(n[2:end], l[2:end], model, f_model[2:end])
    end
  end
  x = Cluster(u, false, convert(CellAdhesionFloat,NaN), f_model[1], n[1], l[1])

  state!(x)

  return x

end


# These functions convert the inputs in the correct data type (CellAdhesionInt and CellAdhesionFloat)

function Cluster(n::Vector{N}, l::Vector{M}, model::T, f_model::Vector{Symbol}) where {T<:BondModel, N<:Real, M<:Real}
  
  n = convert(Vector{CellAdhesionInt}, n)
  l = convert(Vector{CellAdhesionFloat}, l)

  Cluster(n, l, model, f_model)

end

function Cluster(n::N, l::M, model::T, f_model::Symbol) where {T<:BondModel, N<:Real, M<:Real}
  
  n = convert(CellAdhesionInt, n)
  l = convert(CellAdhesionFloat, l)

  Cluster(n, l, model, f_model)

end




function Bond(model::BondModel{T}) where T
  # set state randomly based on the equilibrium binding probability at 0 force.
  v = isless(rand(),k_on(model) / (k_on(model)+k_off(model,0.)))
  return Bond(v, convert(CellAdhesionFloat,NaN), model)
end


# function Bond(model::SlipBondModel)

#   K = model.k_on[:k_on_0] / (model.k_on[:k_on_0] + model.k_off[:k_off_0])
#   v = isless(rand(),K)

#   return Bond(v, convert(CellAdhesionFloat,NaN), model)

# end

# function Bond(model::CatchBondModel)

#   K = model.k_on[:k_on_0] / (model.k_on[:k_on_0] + model.k_off[:k_off_0s]+ model.k_off[:k_off_0c])
#   v = isless(rand(),K)

#   return Bond(v, convert(CellAdhesionFloat,NaN), model)

# end




"""
  bond_state_force(v; output=:flat, time=missing)

Return bond state/force information for a `Bond`, a `Cluster{Bond}`, or nested `Cluster`.

- `output=:flat` (default) returns tuples `(states, forces)`.
- `output=:nested` returns a `Dict(time => force_tree)` representation of the hierarchy,
  where each cluster level stores `[cluster_force, child_forces]`.

Closed bonds are assigned `NaN` in the force output.
"""
function bond_state_force(v::Bond; output::Symbol = :flat, time::Union{Missing,Real} = missing)
  if output == :flat
    state = v.state
    force = v.state ? v.f : convert(CellAdhesionFloat, NaN)
    return Bool[state], CellAdhesionFloat[force]
  elseif output == :nested
    force = v.state ? v.f : convert(CellAdhesionFloat, NaN)
    return Dict(time => force)
  else
    throw(ArgumentError("output must be :flat or :nested"))
  end
end

function bond_state_force(v::Cluster{Bond{T}}; output::Symbol = :flat, time::Union{Missing,Real} = missing) where T <: BondModel
  if output == :flat
    states = getfield.(v.u, :state)
    forces = getfield.(v.u, :f)
    force_out = Vector{CellAdhesionFloat}(undef, v.n)
    for i = 1:v.n
      force_out[i] = states[i] ? forces[i] : convert(CellAdhesionFloat, NaN)
    end
    return states, force_out
  elseif output == :nested
    force_out = Vector{CellAdhesionFloat}(undef, v.n)
    for i = 1:v.n
      force_out[i] = v.u[i].state ? v.u[i].f : convert(CellAdhesionFloat, NaN)
    end
    return Dict(time => Any[v.f, force_out])
  else
    throw(ArgumentError("output must be :flat or :nested"))
  end
end

function bond_state_force(v::Cluster; output::Symbol = :flat, time::Union{Missing,Real} = missing)
  if output == :flat
    states = Bool[]
    forces = CellAdhesionFloat[]
    for i = 1:v.n
      sub_states, sub_forces = bond_state_force(v.u[i]; output = :flat)
      append!(states, sub_states)
      append!(forces, sub_forces)
    end
    return states, forces
  elseif output == :nested
    nested_forces = Vector{Any}(undef, v.n)
    for i = 1:v.n
      child_force_dict = bond_state_force(v.u[i]; output = :nested, time = time)
      nested_forces[i] = child_force_dict[time]
    end
    return Dict(time => Any[v.f, nested_forces])
  else
    throw(ArgumentError("output must be :flat or :nested"))
  end
end




"""
  _save_bond_to_dict(bond::Bond)::Dict

Serialize a Bond's parameters to a dictionary. Includes model type, and all model parameters.
Handles NaN values by converting them to the string "nan" for JSON compatibility.
"""
function _save_bond_to_dict(bond::Bond)::Dict
  # Extract the BondType from the BondModel{T} type parameter
  bond_model_type = typeof(bond.model)
  model_type = String(nameof(bond_model_type.parameters[1]))
  
  # Extract model parameters from the NamedTuple using getfield to access the actual field
  params_tuple = getfield(bond.model, :p)
  model_params = Dict(String(k) => v for (k, v) in pairs(params_tuple))
    
  return Dict(
    "type" => "bond",
    "model_type" => model_type,
    "model_params" => model_params
  )
end


"""
  _save_cluster_to_dict(cluster::Cluster)::Dict

Recursively serialize a Cluster's parameters to a dictionary. Includes cluster-level parameters
(f_model, n, l) and recursively saves all children (which can be Bonds or nested Clusters).
"""
function _save_cluster_to_dict(cluster::Cluster)::Dict  
  # Recursively save children
  children = []
  for i = 1:cluster.n
    child = cluster.u[i]
    if isa(child, Bond)
      push!(children, _save_bond_to_dict(child))
    else  # isa(child, Cluster)
      push!(children, _save_cluster_to_dict(child))
    end
  end
  
  return Dict(
    "type" => "cluster",
    "f_model" => String(cluster.f_model),
    "n" => Int(cluster.n),
    "l" => Float64(cluster.l),
    "children" => children
  )
end


"""
  save_cluster_to_json(cluster::Cluster, filename::String)

Save all parameters of a cluster (or cluster of clusters) to a JSON file.

Example:
```julia
model = SlipBondModel((k_on_0=1.0,), (k_off_0=2.0, f_1e=3.0))
cluster = Cluster(5, 1.0, model, :force_global)
save_cluster_to_json(cluster, "cluster_params.json")
```
"""
function save_cluster_to_json(cluster::Cluster, filename::String)
  dict = Dict()

  # save cluster structure information: cluster has f_model, n, l, Bonds have model type and parameters.
  dict["cluster"] = _save_cluster_to_dict(cluster)
    
  open(filename, "w") do io
    JSON.json(io, dict; allownan=true, nan="nan")
  end
end


"""
  _load_bond_from_dict(bond_dict::Union{Dict, AbstractDict})::Bond

Deserialize a Bond from a dictionary. Reconstructs the BondModel with the correct type and parameters.
"""
function _load_bond_from_dict(bond_dict::Union{Dict, AbstractDict})::Bond
  model_type_str = bond_dict["model_type"]
  model_params = bond_dict["model_params"]
  
  # Reconstruct the BondModel with the saved type and parameters
  if model_type_str == "Slip"
    model_type_class = Slip
  elseif model_type_str == "Catch"
    model_type_class = Catch
  else
    throw(ArgumentError("Unknown bond model type: $model_type_str"))
  end
  
  # Convert parameter names from strings to symbols and values to CellAdhesionFloat
  kwargs = Dict(Symbol(k) => convert(CellAdhesionFloat, v) for (k, v) in model_params)
  
  # Create the BondModel using the generic constructor
  model = BondModel{model_type_class}(; kwargs...)
  
  return Bond(false, convert(CellAdhesionFloat, NaN), model)
end


"""
  _load_cluster_from_dict(cluster_dict::Union{Dict, AbstractDict})::Cluster

Recursively deserialize a Cluster from a dictionary. Reconstructs the complete cluster hierarchy
by recursively loading children (which can be Bonds or nested Clusters).
"""
function _load_cluster_from_dict(cluster_dict::Union{Dict, AbstractDict})::Cluster
  f_model = Symbol(cluster_dict["f_model"])
  n = convert(CellAdhesionInt, cluster_dict["n"])
  l = convert(CellAdhesionFloat, cluster_dict["l"])
  
  children_dicts = cluster_dict["children"]
  
  # Recursively load children and create a vector of the appropriate type
  children = []
  for child_dict in children_dicts
    if child_dict["type"] == "bond"
      push!(children, _load_bond_from_dict(child_dict))
    else  # child_dict["type"] == "cluster"
      push!(children, _load_cluster_from_dict(child_dict))
    end
  end
  
  # Create a Cluster with the loaded data
  cluster = Cluster(children, false, convert(CellAdhesionFloat, NaN), f_model, n, l)
  
  return cluster
end


"""
  load_from_json(filename::String, load_force::Bool)::Cluster

Load a cluster (and force) from a JSON file that was saved with save_cluster_to_json/.

Example:
```julia
cluster = load_from_json("cluster_params.json", false)
```
"""
function load_from_json(filename::String, load_force::Bool)
  data = JSON.parsefile(filename; allownan=true, nan="nan")
  cluster = _load_cluster_from_dict(data["cluster"])
  if load_force
    # get overall force
    force_on_junction = convert(Vector{CellAdhesionFloat}, data["force"])
    # get time step
    dt = convert(CellAdhesionFloat, data["dt"])

    # JSON object keys are strings; keep only numeric time keys and map them to Float64.
    force_time_dict = Dict{Float64, Any}()
    for key in keys(data)
      key_str = String(key)
      t = tryparse(Float64, key_str)
      if !isnothing(t)
        force_time_dict[t] = data[key]
      end
    end

    return cluster, force_on_junction, dt, force_time_dict
  end
  return cluster
end


