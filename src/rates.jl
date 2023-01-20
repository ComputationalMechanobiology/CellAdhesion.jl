export k_off_slip, k_on_constant

"""
k_off_slip(junction::Interface, model::Model)

  Compute the probability of unbinding of each closed bonds within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface data containing the state and parameters of the junction
    - model: Model varible containing the unbinding parameters (k_off_0 = rate of unbinding with no force applied, 
                                                                f_1e = threshold force exponential decay)
  Output parameters:
    - junction: Update the probability of unbinding of each bond
"""

function k_off_slip(junction::Interface, model::Model)

  k_off_0 = model.k_off["k_off_0"]
  f_1e = model.k_off["f_1e"]
  k_off = k_off_0 .* exp.(junction.f ./ f_1e) .* junction.v;

  return Interface(junction.state, junction.n, junction.v, junction.k_on, k_off, junction.f, junction.history)

end


"""
k_on_constant(junction::Interface, model::Model)

  Compute the probability of binding of each open bond within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface data containing the state and parameters of the junction 
    - model: Model varible containing the binding parameters (k_on_0 = rate of binding with no force applied)
  Output parameters:
    - junction: Update the probability of binding of each bond
"""

function k_on_constant(junction::Interface, model::Model)

  k_on_0 = model.k_on["k_on_0"]
  k_on = k_on_0 *(ones(junction.n) .- junction.v);

  return Interface(junction.state, junction.n, junction.v, k_on, junction.k_off, junction.f, junction.history)

end
