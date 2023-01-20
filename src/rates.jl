export k_off_slip, k_on_constant

"""
k_off_slip(junction::Interface, param::Dict)

  Compute the probability of unbinding of each closed bonds within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface data containing the state and parameters of the junction
  Output parameters:
    - junction: Update the probability of unbinding of each bond
"""

function k_off_slip(junction::Interface, param::Dict)

  k_off_0 = param["k_off_0"]
  f_1e = param["f_1e"]
  k_off = k_off_0 .* exp.(junction.f ./ f_1e) .* junction.v;

  return Interface(junction.state, junction.n, junction.v, junction.k_on, k_off, junction.f, junction.history)

end


"""
k_on_constant(junction::Interface, param::Dict)

  Compute the probability of binding of each open bond within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface data containing the state and parameters of the junction 
  Output parameters:
    - junction: Update the probability of binding of each bond
"""

function k_on_constant(junction::Interface, param::Dict)

  k_on_0 = param["k_on_0"]
  k_on = k_on_0 *(ones(junction.n) .- junction.v);

  return Interface(junction.state, junction.n, junction.v, k_on, junction.k_off, junction.f, junction.history)

end
