
export print_cluster, plot_cluster


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
function print_cluster(v::Bond) 

  print("state = ", v.state, ", force = ", v.f, "\n")
  print("model = ", v.model, "\n")
  print("--- \n")

end

function print_cluster(v::Cluster)
  print("********** \n")
  print("Cluster type: ", typeof(v), "\n")
  print("Type of Units = ", typeof(v.u), "\n")
  print("State = ", v.state, ", force = ", v.f, "\n")
  print("Force model = ", v.f_model, ", n = ", v.n, ", l = ", v.l, "\n")
  print("********** \n")

  for i = 1:1:v.n

    k = v.u[i]
    print_cluster(k)

    
  end


end  






# optimised plot_cluster function
function plot_cluster(v::Cluster{Bond{T}}, p, x, y) where T <: BondModel
  states = [bond.state for bond in v.u]
  colors = ifelse.(states, :white, :black)
  y_offsets = y .+ (1:v.n) .* v.l
  scatter!(p, repeat([x], v.n), y_offsets, color = colors, markerstrokecolor = colors, markershape = :square, label="")
end


function plot_cluster(v::Cluster, p, x)
  for i = 1:1:v.n
    k = v.u[i]
    plot_cluster(k, p, x, 1 + (i-1)*v.l)
  end
end 

