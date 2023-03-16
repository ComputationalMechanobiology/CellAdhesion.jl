export interface  #, one_step, junction_simulation


function interface(n::Union{Int, Vector{Int}}, l::Union{Float64, Vector{Float64}}, F::Float64, history::Union{Bool, Vector{Bool}}, model::Model)

    check_length = length(n)
    @assert length(l) == check_length "Length missmatch initialisation"
    @assert length(history) == check_length "Length missmatch initialisation"


    K = model.k_on["k_on_0"] / (model.k_on["k_on_0"] + model.k_off["k_off_0"])
  
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
  
    update_state(x)
    force(x, model)
  
    return x
  
  end
  



# function one_step(junction::Interface, model::Model, s::CellAdhesionFloat)

#     if junction.state == true
#         print("Broken junction")
#         return junction
#     else
#         v = KineticMonteCarlo(junction.v, junction.n, junction.k_on, junction.k_off, model)
#         state = check_state(v)
#         f = force(v, junction.n, model, convert(CellAdhesionFloat, s))

#         if junction.history == false
#             history = false
#         else junction.history != false
#             history = hcat(junction.history, v)          
#         end

#         junction = Interface(state, junction.n, v, junction.k_on, junction.k_off, f, history)

#         junction = k_rate_junction(junction, model, model.k_on["model"])
#         junction = k_rate_junction(junction, model, model.k_off["model"])


#         return junction

#     end

# end




# function junction_simulation(junction::Interface, model::Model, stress::Union{Array{Float64}, Float64}; max_steps::Integer = 1000, verbose::Bool = false)

#     step = 0
#     stress_break = 0
#     time_break = 0

#     # Constant stress
#     if typeof(stress) == Float64

#         while (step <= max_steps) && (junction.state == false)
#             step = step + 1
#             junction = one_step(junction, model, convert(CellAdhesionFloat, stress))
#         end

#         stress_break = stress

#     else 
#         # Arbitrary stress history applied to the junction
#         @assert max_step<=length(stress) "Maximum number of steps exceed stress vector length"

#         while (step <= max_steps) && (junction.state == false)
#             step = step + 1
#             junction = one_step(junction, model, convert(CellAdhesionFloat, stress[step]))
#         end

#         stress_break = stress[step]

#     end

#     if verbose == true
#         if junction.state == true
#             print("Junction broken")
#         elseif step > max_steps
#             print("Maximum number of iterations reached")
#         end
#     end

#     time_break = model.param["dt"]*step

#     return junction, stress_break, time_break, step


# end