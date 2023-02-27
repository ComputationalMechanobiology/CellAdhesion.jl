# export init_junction, one_step, junction_simulation

# function init_junction(n::Integer, model::Model; history::Bool = false)


#     K = model.k_on["k_on_0"] / (model.k_on["k_on_0"] + model.k_off["k_off_0"])

#     # if !haskey(model.param, "load")
#     #     model.param["load"] = "global"
#     # end
#     # if !haskey(model.param, "dt")
#     #     model.param["dt"] = 0.1
#     # end

#     v = init_bonds(n, convert(CellAdhesionFloat, K))
#     state = check_state(v)
#     f = force(v, n, model, convert(CellAdhesionFloat, 0))

#     if history == false
#         junction = Interface(state, n, v, [], [], f, false)
#     else
#         junction = Interface(state, n, v, [], [], f, v)
#     end

#     junction = k_rate_junction(junction, model, model.k_on["model"])
#     junction = k_rate_junction(junction, model, model.k_off["model"])

#     return junction

# end



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