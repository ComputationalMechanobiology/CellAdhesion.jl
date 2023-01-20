
export init_bonds, check_state, force

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
check_state(junction::Interface)

  Check if a junction is broken or still viable.

  Input parameters:
    - junction: Interface data containing the state of each single bond (junction.v)
  Output parameters:
    - junction: Update the state of the junction (true = broken, false = viable)
"""

function check_state(junction::Interface)

  @assert !isempty(junction.v) "Bond state vector in Interface is empty"

  sum_v = sum(junction.v)
  state = isequal(sum_v,0)
  
  return Interface(state, junction.n, junction.v, junction.k_on, junction.k_off, junction.f, junction.history)

end


"""
distance(junction::Interface)

  Compute the distance of each from its two closest closed link on each side. 
  Used to compute a "local" load distributions across closed bonds.
  Periodic bondary conditions for the edges.

  Input parameters:
    - junction: Interface data containing the state and parameters of the junction 
  Output parameters:
    - l: vector of CellAdhesionFloat containing the distance for each closed bond
"""

function distance(junction::Interface)

  @assert !isempty(junction.v) "Bond state vector in Interface is empty"

  l = zeros(junction.n);
  if sum(junction.v) == 0
       l .= zeros(junction.n);
  elseif sum(junction.v) ==1
      temp = findall(x->x==1, junction.v)
      l[temp[1]] = junction.n;        
  elseif sum(junction.v) == 2
      temp = findall(x->x==1, junction.v)
      l[temp[1]] = junction.n;
      l[temp[2]] = junction.n;    
  else
      temp = findall(x->x==1, junction.v)
      gaps = diff(temp);
      dist = zeros(size(temp));
      dist[2:end-1] = gaps[2:end] .+ gaps[1:end-1];
      dist[1] = junction.n - temp[end] + temp[2];
      dist[end] = junction.n - temp[end-1] + temp[1];
      l[temp] .= dist;
  end

  return convert(Vector{CellAdhesionFloat}, l)

end


"""
force(junction::Interface, param::Dict, s::CellAdhesionFloat)

  Compute the force on each link. 
  Two options are available:
   - global: equal distribution across all closed bonds
   - local: inhomogeneous distribution of load accounting fro the distance to its nearest closed neighbor on both sides. 

   If no load parameter is specified, global is used by default

  Input parameters:
    - junction: Interface data containing the state of the junction 
    - param: parameters of the junction (containing the type of load distribution: "local" or "global")
    - s: applied stress to the junction (CellAdhesionFloat type)
  Output parameters:
    - junction: updated force applied to each bond within the Interface data structure
"""

function force(junction::Interface, param::Dict, s::CellAdhesionFloat)

  @assert !isempty(junction.v) "Bond state vector in Interface is empty"

  alpha = zeros(junction.n);
  if haskey(param, "load")
    state = param["load"]
  else
    state = "global"
  end

  if param["load"] == "global"
      alpha .= junction.n/sum(junction.v) .*junction.v;
  elseif param["load"] == "local"
      l = distance(junction)
      alpha .= junction.n .* l ./ sum(l);
  end

  return Interface(junction.state, junction.n, junction.v, junction.k_on, junction.k_off, convert(Vector{CellAdhesionFloat}, alpha .* s), junction.history)

end







function KineticMonteCarlo!(v,k_off, k_on,N,k_off_0, dt)

  k_tot = k_off.+k_on;
  random = rand(N);

  bond_events = findall(x->x==1, ((k_on.*dt) .>random));
  v[bond_events] .= 1;
  unbond_events = findall(x->x==1, ((k_off.*dt) .>random));
  v[unbond_events] .= 0;


end


