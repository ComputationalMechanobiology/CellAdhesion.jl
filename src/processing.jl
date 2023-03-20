export interface, one_step, force_increment, junction_simulation


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
    k_rate_junction(x, model)
  
    return x
  
  end
  



function one_step(v::Interface, model::Model)

  update_state(v)

  if v.state == false
      print("Broken junction")
  else
      KineticMonteCarlo(v,model)
      update_state(v)
      force(v, model)
      k_rate_junction(v, model)
  end

end

function force_increment(v::Interface, model::Model, F::CellAdhesionFloat)

  setfield!(v, :f, F)
  force(v, model)
  k_rate_junction(v, model)


end




function junction_simulation(junction::Interface, model::Model, force::Union{Array{Float64}, Float64}; max_steps::Integer = 1000, verbose::Bool = false)

    step = 0
    force_break = 0
    time_break = 0

    if typeof(force) == Float64
      force = convert(CellAdhesionFloat,force)
    else
      force = convert(Vector{CellAdhesionFloat},force)
    end

    # Constant force
    if typeof(force) == CellAdhesionFloat
      force_increment(junction, model, force)

      #print(junction.u, "\n")

      while (step <= max_steps) && (junction.state == true)
          step = step + 1
          one_step(junction, model)
      end

      force_break = force

    else 

      # Arbitrary force history applied to the junction
      @assert max_step<=length(force) "Maximum number of steps exceed force vector length"

      while (step <= max_steps) && (junction.state == true)
          step = step + 1
          force_increment(junction, model, force[i])
          one_step(junction, model)
      end

      force_break = force[step]

    end

    if verbose == true
        if junction.state == false
            print("Junction broken")
        elseif step > max_steps
            print("Maximum number of iterations reached")
        end
    end

    time_break = model.param["dt"]*step


    return force_break, time_break, step


end