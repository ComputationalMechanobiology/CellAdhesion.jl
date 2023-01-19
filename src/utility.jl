
export init_bonds, check_state, link_k_off_slip, link_k_on_constant

"""
init_bond(n::Integer, K::CellAdhesionFloat)

  Generate the initial state of each bond within the cell junction.

  Input parameters:
    - n: number of bonds in the junction
    - K: probability to start in a closed state
  Output parameters:
    - v: state vecotr of individual bonds (0 = open, 1 = closed)
"""

function init_bonds(n::Integer,K::CellAdhesionFloat)

    @assert n>0 "Number of bonds cannot be negative or equal to 0"
    @assert (K>=0) && (K<=1) "Initialisation probability must be 0<=K<=1"

    v = isless.(rand(n),K)


    return v

end


"""
check_state(v)

  Check if a junction is broken or still viable.

  Input parameters:
    - junction: Interface data containing the state of each single bond (junction.v)
  Output parameters:
    - junction: Update the state of the junction (true = broken, false = viable)
"""

function check_state(junction::Interface)


  sum_v = sum(junction.v)
  state = isequal(sum_v,0)
  
  return Interface(state, junction.n, junction.v, junction.k_on, junction.k_off, junction.f, junction.param)

end


"""
link_k_off_slip(v)

  Compute the probability of unbinding of each closed bonds within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface data containing the state and parameters of the junction
  Output parameters:
    - junction: Update the probability of unbinding of each bond
"""

function link_k_off_slip(junction::Interface)

  k_off_0 = junction.param["k_off_0"]
  f_1e = junction.param["f_1e"]
  k_off = k_off_0 .* exp.(junction.f ./ f_1e) .* junction.v;

  return Interface(junction.state, junction.n, junction.v, junction.k_on, k_off, junction.f, junction.param)

end


"""
link_k_on_constant(v)

  Compute the probability of binding of each open bond within the junction. 
  Force-dependent behaviour of the unbinding probability is described by the Bell model. 

  Input parameters:
    - junction: Interface data containing the state and parameters of the junction 
  Output parameters:
    - junction: Update the probability of binding of each bond
"""

function link_k_on_constant(junction::Interface)

  k_on_0 = junction.param["k_on_0"]
  k_on = k_on_0 *(ones(junction.n) .- junction.v);

  return Interface(junction.state, junction.n, junction.v, k_on, junction.k_off, junction.f, junction.param)

end