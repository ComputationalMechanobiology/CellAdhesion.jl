export setforce!

function k_on(m::SlipBondModel)
  
  return m.k_on[:k_on_0]

end

function k_off(m::SlipBondModel, f::CellAdhesionFloat)
  
  return m.k_off[:k_off_0] .* exp.(f ./ m.k_off[:f_1e])

end

#------------------ FORCE -------------------------------------

"""
setforce!(v::Cluster{Bond{T}}, F::CellAdhesionFloat)

  Input parameters:
    - v: Interface structure
    - F: Force applied to the cluster
  Output parameters:
    - Updated Interface with force applied to each link
"""
function setforce!(v::Cluster{Bond{T}}, F::CellAdhesionFloat) where T <:BondModel

  setfield!(v, :f, F)
  distributeforce!(v)

end


function setforce!(v::Cluster{Bond{T}}) where T <:BondModel
  
  distributeforce!(v)

end


#function setforce!(v::Cluster, F::CellAdhesionFloat)

#  setfield!(v, :f, F)
#  setforce!(v)

#end


function setforce!(v::Cluster)

  if v.state == true

    distributeforce!(v)

    for i = 1:1:v.n 
      k = v.u[i]
      setforce!(k)  
    end

  end

end


function distributeforce!(v::Cluster)

  #@assert v.f>=0 "Applied stress to Cluster must be positive or equal to zero"
  update_f = getfield(CellAdhesion, v.f_model)(v)      # If we define v.f_model as string then => Symbol(v.f_model)
  for i = 1:1:v.n
    setfield!(v.u[i], :f, update_f[i])
  end

  return update_f
end


"""
force_global
Computer force distribution by equally dividing the force within the closed bonds
"""
function force_global(v::Cluster)

  interface_v = getfield.(v.u, :state);
  
  return interface_v .* v.f./sum(interface_v)
  
end

"""
force_local
Computer force distribution by accounting for the distance of each link from its two closest neighbours
"""
function force_local(v::Cluster)

  interface_v = getfield.(v.u, :state);
  l = distance(interface_v, v.n)
  return l ./ sum(l) .*v.f;

end

"""
distance(v::BitVector, n::Integer)

  Compute the distance of each from its two closest closed link on each side. 
  Used to compute a "local" load distributions across closed bonds.
  Periodic bondary conditions for the edges.

  Input parameters:
    - v: vector with the state of each single bond
    - n: number of bonds in the junction
  Output parameters:
    - l: vector of CellAdhesionFloat containing the distance for each closed bond
"""
function distance(v::BitVector, n::Integer)

  @assert !isempty(v) "Bond state vector in Interface is empty"

  l = zeros(n);
  if sum(v) == 0
       l .= zeros(n);
  elseif sum(v) ==1
      temp = findall(x->x==1, v)
      l[temp[1]] = n;        
  elseif sum(v) == 2
      temp = findall(x->x==1, v)
      l[temp[1]] = n;
      l[temp[2]] = n;    
  else
      temp = findall(x->x==1, v)
      gaps = diff(temp);
      dist = zeros(size(temp));
      dist[2:end-1] = gaps[2:end] .+ gaps[1:end-1];
      dist[1] = n - temp[end] + temp[2];
      dist[end] = n - temp[end-1] + temp[1];
      l[temp] .= dist;
  end

  return convert(Vector{CellAdhesionFloat}, l)

end
