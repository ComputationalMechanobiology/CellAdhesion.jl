export update!, runcluster, Cluster, Bond, bond_state_force, save_cluster_state


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

  #if v.state == true

    for i = 1:1:v.n 
      k = v.u[i]
      update!(k, dt)  
    end

    #Get the state value for each bond
    interface_v = getfield.(v.u, :state);

    # If the sum of the state values is 0, the junction is broken 
    sum_v = sum(interface_v);
    state = isequal(sum_v,0);

    # Update the state value of the junction
    setfield!(v, :state, !state)

  # end

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

function runcluster(v::Cluster, force::Vector{Float64}, dt::Float64, file_name::String, ; max_steps::Integer = 1000, verbose::Bool = false)
  # dispatch to the appropriate method based on file extension
  ext = splitext(file_name)[2]
  if ext == ".json"
    return runcluster_json(v, force, dt, file_name; max_steps = max_steps, verbose = verbose)
  elseif ext in [".svg", ".png", ".jpg", ".jpeg"]
    return runcluster_plot(v, force, dt, file_name; max_steps = max_steps, verbose = verbose)
  else
    throw(ArgumentError("Unsupported file extension: $ext. Supported extensions are .json, .svg, .png, .jpg, .jpeg"))
  end
end 


function runcluster_plot(v::Cluster, force::Vector{Float64}, dt::Float64, plot_file_name::String, ; max_steps::Integer = 1000, verbose::Bool = false)

  # Arbitrary force history applied to the junction
  n = length(force)

  if max_steps > n
         @warn max_steps<=n "Maximum number of steps exceed force vector length"
  	 max_steps = n
         print("\n Maximum number of steps = ", max_steps, "\n")
  end
  

  step = 1
  p = plot()
  plot_cluster(v, p, step, 0)

  force = convert(Vector{CellAdhesionFloat},force)
  dt = convert(CellAdhesionFloat, dt)

  while (step <= max_steps) && (v.state == true)
      F = force[step]
      setforce!(v, F)
      update!(v, dt)
      step = step + 1
      plot_cluster(v, p, step, 0)
  end

  if verbose == true
      if v.state == false
          print("Junction broken")
      elseif step > max_steps
          print("Maximum number of iterations reached")
      end
  end

  savefig(p, plot_file_name)

  return v.state, force[step-1], dt*(step-1), (step-1)
end



function runcluster_json(v::Cluster, force::Vector{Float64}, dt::Float64, json_file_name::String; max_steps::Integer = 1000, verbose::Bool = false)
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
  save_cluster_state(states, json_file_name)
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
  x = Cluster(u, false, convert(CellAdhesionFloat, 0.0), f_model, n, l)
  
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
  x = Cluster(u, false, convert(CellAdhesionFloat, 0.0), f_model[1], n[1], l[1])

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



function Bond(model::SlipBondModel)

  K = model.k_on[:k_on_0] / (model.k_on[:k_on_0] + model.k_off[:k_off_0])
  v = isless(rand(),K)

  return Bond(v, convert(CellAdhesionFloat, 0.0), model)

end

function Bond(model::CatchBondModel)

  K = model.k_on[:k_on_0] / (model.k_on[:k_on_0] + model.k_off[:k_off_0s]+ model.k_off[:k_off_0c])
  v = isless(rand(),K)

  return Bond(v, convert(CellAdhesionFloat, 0.0), model)

end




"""
  bond_state_force(v; output=:flat, time=missing)

  Return bond state/force information for a `Bond`, a `Cluster{Bond}`, or nested `Cluster`.

  - `output=:flat` (default) returns tuples `(states, forces)` identical to the previous API.
  - `output=:nested` returns a `Dict` with keys "time", "states", "force" that includes
    each level of the hierarchy. At every cluster level, the first entry is the cluster
    state/force, followed by its children. Example structure:

    * states = [true, [[true, [true, true, false]],
                       [false, [false, false, false]],
                       [true, [true, false, false]]]]
    * force  = [4.0, [[3.0, [1.0, 2.0, NaN]],
                      [NaN, [NaN, NaN, NaN]],
                      [1.0, [1.0, NaN, NaN]]]]

  Closed bonds are assigned `NaN` in the force output.
"""

function _state_force_tree(v::Bond)
  force = v.state ? v.f : convert(CellAdhesionFloat, NaN)
  return v.state, force
end

function _state_force_tree(v::Cluster{Bond{T}}) where T <: BondModel
  states = getfield.(v.u, :state)
  force_out = Vector{CellAdhesionFloat}(undef, v.n)
  for i = 1:v.n
    force_out[i] = v.u[i].state ? v.u[i].f : convert(CellAdhesionFloat, NaN)
  end
  return Any[v.state, states], Any[v.f, force_out]
end

function _state_force_tree(v::Cluster)
  nested_states = Vector{Any}(undef, v.n)
  nested_forces = Vector{Any}(undef, v.n)
  for i = 1:v.n
    nested_states[i], nested_forces[i] = _state_force_tree(v.u[i])
  end
  return Any[v.state, nested_states], Any[v.f, nested_forces]
end

function _state_force_dict(v, time)
  states_tree, force_tree = _state_force_tree(v)
  return Dict(
    "time" => time,
    "states" => states_tree,
    "force" => force_tree,
  )
end

function bond_state_force(v::Bond; output::Symbol = :flat, time::Union{Missing,Real} = missing)
  if output == :flat
    state = v.state
    force = v.state ? v.f : convert(CellAdhesionFloat, NaN)
    return Bool[state], CellAdhesionFloat[force]
  elseif output == :nested
    return _state_force_dict(v, time)
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
    return _state_force_dict(v, time)
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
    return _state_force_dict(v, time)
  else
    throw(ArgumentError("output must be :flat or :nested"))
  end
end



"""
  save_cluster_state(v, file_name; output=:nested, time=missing)
  save_cluster_state(states, file_name)

Save cluster bond state/force information to a JSON file.

First method:
- Input: `v::Cluster`, `file_name::String`.
- Optional kwargs: `output::Symbol` (default `:nested`), `time::Union{Missing,Real}` (default `missing`).
- Uses `bond_state_force(v; output=output, time=time)` to build a state dictionary, then writes it as JSON.

Second method:
- Input: `states::Vector{<:AbstractDict}`, `file_name::String`.
- Writes a vector of already-obtained bond state dictionaries to JSON.

Also supported helper overload:
- `save_cluster_state(state::Dict, file_name::String)`.

JSON encoding handles nested arrays/dicts, booleans, numbers, strings, and missing values.
"""

function _save_json_to_file(obj, file_name::String)
  open(file_name, "w") do io
    JSON.json(io, obj; allownan=true, nan="nan")
  end
  return file_name
end

function save_cluster_state(v::Cluster, file_name::String; output::Symbol = :nested, time::Union{Missing,Real} = missing)
    state = bond_state_force(v; output = output, time = time)
    _save_json_to_file(state, file_name)
end

function save_cluster_state(states::Vector{<:AbstractDict}, file_name::String)
    _save_json_to_file(states, file_name)
end

# convenience additional overload
function save_cluster_state(state::Dict, file_name::String)
    _save_json_to_file(state, file_name)
end


