export k_off_slip, k_on_constant



function k_rate_junction(junction, model, k_function)

    @assert (k_function in list_k_models) "k rate model not defined"
    junction = k_function(junction, model)  
    return junction

end


"""
k_off_slip(junction::Interface, model::Model)

  Compute the probability of unbinding of each closed bonds within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface variable
    - model: Model varible containing the unbinding parameters (k_off_0 = rate of unbinding with no force applied, 
                                                                f_1e = threshold force exponential decay)
  Output parameters:
    - k_off: Unbinding probability for each bond
"""

function k_off_slip(junction::Interface, model::Model)

   k_off = model.k_off["k_off_0"] .* exp.(junction.f ./ model.k_off["f_1e"]) .* junction.v;

   return Interface(junction.state, junction.n, junction.v, junction.k_on, convert(Vector{CellAdhesionFloat}, k_off), junction.f, junction.history)

end


"""
k_on_constant(junction::Interface, model::Model)

  Compute the probability of binding of each open bond within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface variable
    - model: Model varible containing the binding parameters (k_on_0 = rate of binding with no force applied)
  Output parameters:
    - junction: binding probability for each bond
"""

function k_on_constant(junction::Interface, model::Model)

  k_on = model.k_on["k_on_0"] *(ones(junction.n) .- junction.v);
  
  return Interface(junction.state, junction.n, junction.v, convert(Vector{CellAdhesionFloat}, k_on), junction.k_off, junction.f, junction.history)

end

list_k_models = [k_off_slip, 
                 k_on_constant]