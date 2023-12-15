
export print_cluster


"""
state!(v)

Updates the state field of a Cluster or a Bond.
When all the subunits are open, then the state of the structure is set to open. Otherwise, it is close.   

Input paramters:
  - v: interface that can either be a Cluster or a Bond. 

"""
function state!(v::Bond)

  return getfield(v, :state)

end

function state!(v::Cluster)

  interface_v = Vector{Bool}(undef,v.n)

  for i = 1:1:v.n 
    k = v.u[i]
    interface_v[i] = state!(k)  
  end

  # If the sum of the state values is 0, the junction is broken 
  sum_v = sum(interface_v);
  state = isequal(sum_v,0);

  # Update the state value of the junction
  setfield!(v, :state, !state)
  

end



"
Base.setproperty!   

Definition to update Interface struct fields
"

function Base.setproperty!(x::Cluster, s::Symbol, new_x::Vector{CellAdhesionFloat})

  for i = 1:1:x.n
    setfield!(x.u[i], s, new_x[i])
  end

end




"""
  print_cluster(x)


  Nice screen print of Cluster structure or Bond (x)

"""
function print_cluster(x::Bond) 

  print("state = ", x.state, ", force = ", x.f, "\n")
  print("model = ", x.model, "\n")
  print("--- \n")

end

function print_cluster(x::Cluster)
  print("********** \n")
  print("Cluster type: ", typeof(x), "\n")
  print("Type of Units = ", typeof(x.u), "\n")
  print("State = ", x.state, ", force = ", x.f, "\n")
  print("Force model = ", x.f_model, ", n = ", x.n, ", l = ", x.l, "\n")
  print("********** \n")

  for i = 1:1:x.n

    k = x.u[i]
    print_cluster(k)

    
  end


end  
