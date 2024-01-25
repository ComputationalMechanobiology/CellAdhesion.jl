export update!, runcluster, Cluster, Bond


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



function Bond(model::T) where T <: BondModel

  K = model.k_on[:k_on_0] / (model.k_on[:k_on_0] + model.k_off[:k_off_0])
  v = isless(rand(),K)

  return Bond(v, convert(CellAdhesionFloat, 0.0), model)

end

