export update!, cluster_simulation


function update!(v::Bond, dt::CellAdhesionFloat)

  random = rand();
  v.state ? k = k_off(v.model, v.f) : k = k_on(v.model)
  (k*dt>random) ? (v.state = !(v.state)) : nothing

  
end

function update!(v::Cluster, dt::CellAdhesionFloat)

  if v.state == true

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

  end


end



function cluster_simulation(v::Cluster, force::Float64, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false)

  step = 0

  force = convert(CellAdhesionFloat,force)
  dt = convert(CellAdhesionFloat, dt)

  while (step <= max_steps) && (v.state == true)
    step = step + 1
    setforce!(v, force)
    update!(v, dt)
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


function cluster_simulation(v::Cluster, force::Vector{Float64}, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false)

  # Arbitrary force history applied to the junction
  @assert max_step<=length(force) "Maximum number of steps exceed force vector length"

  step = 1
  F = 0

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

  return v.state, F, dt*step, step


end




# function Cluster(n::Int64, l::Float64, model::T, f_model::Symbol) where T<:BondModel

#   @assert n>0 "Number of bonds cannot be negative or equal to 0"
#   n = convert(CellAdhesionInt, n)
#   l = convert(CellAdhesionFloat, l)
  
#   x = init_bonds(n, l, model, F, f_model)

#   return x
    
# end




function Cluster(n::CellAdhesionInt, l::CellAdhesionFloat, model::T, f_model::Symbol) where T<:BondModel

  u = Vector{Bond}(undef, n)

  for i = 1:1:n
    u[i] = Bond(model)
  end
  x = Cluster(u, false, convert(CellAdhesionFloat, 0.0), f_model, n, l)
  check_state!(x)

  return x

end

function Bond(model::T) where T <: BondModel

  K = model.k_on[:k_on_0] / (model.k_on[:k_on_0] + model.k_off[:k_off_0])
  v = isless(rand(),K)

  return Bond(v, convert(CellAdhesionFloat, 0.0), model)

end





# export interface, one_step, force_increment, junction_simulation


# function interface(n::Union{Int, Vector{Int}}, l::Union{Float64, Vector{Float64}}, model::T, F::Union{Int64,Float64}, f_model::Union{Symbol, Vector{Symbol}})

#     check_length = length(n)
#     @assert length(l) == check_length "Length missmatch initialisation"

  
#     # Make variable types consistent with CellAdhesionFloat
#     if check_length == 1
      
#       @assert n>0 "Number of bonds cannot be negative or equal to 0"

#       n = convert(CellAdhesionInt, n)
#       l = convert(CellAdhesionFloat, l)
#     else
#       n = convert(Vector{CellAdhesionInt}, n)
#       l = convert(Vector{CellAdhesionFloat}, l)
#     end

#     F = convert(CellAdhesionFloat, F)

#     x = init_bonds(n, l, model, F, f_model)
  
#     update_state(x)
  
#     return x
  
#   end
  



# function one_step(v::Union{Cluster,Interface}, dt::CellAdhesionFloat)

#   update_state(v)

#   if v.state == false
#       print("Broken junction")
#   else
#       KineticMonteCarlo(v,dt)
#       update_state(v)
#       force(v)
#       k_rate_junction(v)
#   end

# end

# function force_increment(v::Union{Cluster, Interface}, F::CellAdhesionFloat)

#   setfield!(v, :f, F)
#   force(v)
#   k_rate_junction(v)

# end











# function junction_simulation(junction::Union{Cluster,Interface}, force::Union{Array{Float64}, Float64}, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false)

#     step = 0
#     force_break = 0
#     time_break = 0

#     if typeof(force) == Float64
#       force = convert(CellAdhesionFloat,force)
#     else
#       force = convert(Vector{CellAdhesionFloat},force)
#     end

#     dt = convert(CellAdhesionFloat, dt)

#     # Constant force
#     if typeof(force) == CellAdhesionFloat
#       force_increment(junction, force)

#       #print(junction.u, "\n")

#       while (step <= max_steps) && (junction.state == true)
#           step = step + 1
#           one_step(junction, dt)
#       end

#       force_break = force

#     else 

#       # Arbitrary force history applied to the junction
#       @assert max_step<=length(force) "Maximum number of steps exceed force vector length"

#       while (step <= max_steps) && (junction.state == true)
#           step = step + 1
#           force_increment(junction, force[i])
#           one_step(junction, dt)
#       end

#       force_break = force[step]

#     end

#     if verbose == true
#         if junction.state == false
#             print("Junction broken")
#         elseif step > max_steps
#             print("Maximum number of iterations reached")
#         end
#     end

#     time_break = dt*step


#     return force_break, time_break, step


# end