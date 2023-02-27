
export init_bonds, check_state, force


"""
check_state(v::Interface)

  Check if a junction is broken or still viable.

  Input parameters:
    - v: vector with the state of each single bond 
  Output parameters:
    - state: true = broken, false = viable
"""

function check_state(v::Interface)

  @assert !isempty(v.bonds) "Bond vector in Interface is empty"

  interface_v = getfield.(v.bonds, :state)

  sum_v = sum(interface_v)
  state = isequal(sum_v,0)

  setfield!(v, :state, !state)
  

end

























# """
# init_bonds(n::Integer, K::CellAdhesionFloat)

#   Generate the initial state of each bond within the cell junction.

#   Input parameters:
#     - n: number of bonds in the junction
#     - K: probability to start in a closed state
#   Output parameters:
#     - v: state vecotr of individual bonds (0 = open, 1 = closed)
# """

# function init_bonds(n::Integer,K::CellAdhesionFloat)

#     @assert n>0 "Number of bonds cannot be negative or equal to 0"
#     @assert (K>=0) && (K<=1) "Initialisation probability must be 0<=K<=1"

#     v = isless.(rand(n),K)
#     return v

# end



# """
# distance(v::BitVector, n::Integer)

#   Compute the distance of each from its two closest closed link on each side. 
#   Used to compute a "local" load distributions across closed bonds.
#   Periodic bondary conditions for the edges.

#   Input parameters:
#     - v: vector with the state of each single bond
#     - n: number of bonds in the junction
#   Output parameters:
#     - l: vector of CellAdhesionFloat containing the distance for each closed bond
# """

# function distance(v::BitVector, n::Integer)

#   @assert !isempty(v) "Bond state vector in Interface is empty"

#   l = zeros(n);
#   if sum(v) == 0
#        l .= zeros(n);
#   elseif sum(v) ==1
#       temp = findall(x->x==1, v)
#       l[temp[1]] = n;        
#   elseif sum(v) == 2
#       temp = findall(x->x==1, v)
#       l[temp[1]] = n;
#       l[temp[2]] = n;    
#   else
#       temp = findall(x->x==1, v)
#       gaps = diff(temp);
#       dist = zeros(size(temp));
#       dist[2:end-1] = gaps[2:end] .+ gaps[1:end-1];
#       dist[1] = n - temp[end] + temp[2];
#       dist[end] = n - temp[end-1] + temp[1];
#       l[temp] .= dist;
#   end

#   return convert(Vector{CellAdhesionFloat}, l)

# end


# """
# force(v::BitVector, n::Integer, model::Model, s::CellAdhesionFloat)

#   Compute the force on each link. 
#   Two options are available:
#    - global: equal distribution across all closed bonds
#    - local: inhomogeneous distribution of load accounting fro the distance to its nearest closed neighbor on both sides. 

#    If no load parameter is specified, global is used by default

#   Input parameters:
#     - v: vector with the state of each single bond
#     - n: number of bonds in the junction
#     - model: Model structure containing the type of load distribution: "local" or "global" (Model.param)
#     - s: applied stress to the junction (CellAdhesionFloat type)
#   Output parameters:
#     - f: (vector) force applied to each bond 
# """

# function force(v::BitVector, n::Integer, model::Model, s::CellAdhesionFloat)

#   @assert !isempty(v) "Bond state vector in Interface is empty"

#   alpha = zeros(n);
  
#   state = model.param["load"]


#   if state == "global"
#       alpha .= n/sum(v) .*v;
#   elseif state == "local"
#       l = distance(v, n)
#       alpha .= n .* l ./ sum(l);
#   end

#   return convert(Vector{CellAdhesionFloat}, alpha .* s)

# end


# """
# KineticMonteCarlo(v::BitVector, n::Integer, k_on::Vector{CellAdhesionFloat}, k_off::Vector{CellAdhesionFloat}, model::Model)-NOT TESTED YET!!!!!

# """

# function KineticMonteCarlo(v::BitVector, n::Integer, k_on::Vector{CellAdhesionFloat}, k_off::Vector{CellAdhesionFloat}, model::Model)

#   random = rand(n);
#   v_temp = copy(v);

#   bond_events = findall(x->x==1, ((k_on.*model.param["dt"]) .>random));
#   v_temp[bond_events] .= 1;
#   unbond_events = findall(x->x==1, ((k_off.*model.param["dt"]) .>random));
#   v_temp[unbond_events] .= 0;


#   return v_temp


# end


