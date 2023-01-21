export init_junction, one_step

function init_junction(n::Integer, model::Model; history::Bool = false)


    K = model.k_on["k_on_0"] / (model.k_on["k_on_0"] + model.k_off["k_off_0"])

    # if !haskey(model.param, "load")
    #     model.param["load"] = "global"
    # end
    # if !haskey(model.param, "dt")
    #     model.param["dt"] = 0.1
    # end

    v = init_bonds(n, convert(CellAdhesionFloat, K))
    state = check_state(v)
    f = force(v, n, model, convert(CellAdhesionFloat, 0))

    if history == false
        junction = Interface(state, n, v, [], [], f, false)
    else
        junction = Interface(state, n, v, [], [], f, v)
    end

    junction = k_rate_junction(junction, model, model.k_on["model"])
    junction = k_rate_junction(junction, model, model.k_off["model"])

    return junction

end



function one_step(junction::Interface, model::Model, s::Float64)

    if junction.state == true
        print("Broken junction")
        return junction
    else
        v = KineticMonteCarlo(junction.v, junction.n, junction.k_on, junction.k_off, model)
        state = check_state(v)
        f = force(v, junction.n, model, convert(CellAdhesionFloat, s))

        if junction.history == false
            history = false
        else junction.history != false
            history = hcat(junction.history, v)          
        end

        junction = Interface(state, junction.n, v, junction.k_on, junction.k_off, f, history)

        junction = k_rate_junction(junction, model, model.k_on["model"])
        junction = k_rate_junction(junction, model, model.k_off["model"])

        return junction

    end

end




function junction_simulation(junction::Interface, model::Model, stress::Array{Float64}; max_steps::Integer = 1000)











end