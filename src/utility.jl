
# export model_init, update_state, KineticMonteCarlo_unit, KineticMonteCarlo

export slip_model_init, print_cluster


function slip_model_init(k_on::NamedTuple, k_off::NamedTuple)

  k_on_0 = convert(CellAdhesionFloat, k_on[:k_on_0])
  k_off_0 = convert(CellAdhesionFloat, k_off[:k_off_0])
  f_1e = convert(CellAdhesionFloat, k_off[:f_1e])

  return SlipBondModel((k_on_0 = k_on_0,), (k_off_0 = k_off_0, f_1e = f_1e))

end





function check_state!(v::Bond)

  return getfield(v, :state)

end

function check_state!(v::Cluster)

  interface_v = Vector{Bool}(undef,v.n)

  for i = 1:1:v.n 
    k = v.u[i]
    interface_v[i] = check_state!(k)  
  end

  # If the sum of the state values is 0, the junction is broken 
  sum_v = sum(interface_v);
  state = isequal(sum_v,0);

  # Update the state value of the junction
  setfield!(v, :state, !state)
  

end






"""
Base.setproperty!   

Definition to update Interface struct fields
"""


function Base.setproperty!(x::Cluster, s::Symbol, new_x::Vector{CellAdhesionFloat})

  for i = 1:1:x.n
    setfield!(x.u[i], s, new_x[i])
  end

end


"""
print_cluster()

Nice screen print of Cluster structure

"""
# function print_cluster(k::Bond)
#   print("state = ", k.state, ", force = ", k.f, "\n")
#   print("model = ", k.model, "\n")

# end


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



# init_bonds(n::Integer, l::CellAdhesionFloat, F::CellAdhesionFloat, model::Base.RefValue{Model})

#   Generate the initial state of each bond within the cell junction.

#   Input parameters:
#     - n: number of bonds in the junction
#     - 
#   Output parameters:
#     - Interface struct
# """




# function init_bonds(n::Vector{CellAdhesionInt}, l::Vector{CellAdhesionFloat}, model::Base.RefValue{Model}, F::CellAdhesionFloat, f_model::Vector{Symbol})

#   v = Interface(Vector{Cluster}(undef,n[1]), false, F, f_model[1], n[1], l[1])

#   for i = 1:1:n[1]
#     v.u[i] = init_bonds(n[2], l[2], model, F, f_model[2])
#   end

#   return v

# end



# """
# update_unit_state(v::Union{Cluster,Interface})

#   update a cluster or interface state

# """

# function update_unit_state(v::Union{Cluster,Interface})

#   @assert !isempty(v.u) "Unit vector in Cluster is empty"

#   # Get the state value for each bond
#   interface_v = getfield.(v.u, :state);

#   # If the sum of the state values is 0, the junction is broken 
#   sum_v = sum(interface_v);
#   state = isequal(sum_v,0);

#   # Update the state value of the junction
#   setfield!(v, :state, !state)

#   return !state
  

# end


# """
# update_state(v::Interface)
# update_state(v::Cluster)

#   update if an interface or cluster is broken or still viable.

#   Input parameters:
#     - v: Interface object 
#   Output parameters:
#     - v with updated state field
# """


# function update_state(v::Cluster)

#   update_unit_state(v)

# end

# function update_state(v::Interface)

#   for i = 1:1:v.n
#     update_unit_state(v.u[i])
#   end

#   update_unit_state(v)

# end




# function KineticMonteCarlo(v::Interface, dt::CellAdhesionFloat)

#   for i=1:1:v.n
#     KineticMonteCarlo(v.u[i], dt)
#   end

# end




# """
# KineticMonteCarlo(v::Interface, dt::CellAdhesionFloat)

# """

# function KineticMonteCarlo(v::Cluster, dt::CellAdhesionFloat)

#   random = rand(v.n);
#   v_temp = getfield.(v.u, :state);

#   bond_events = findall(x->x==1, (getfield.(v.u, :k_on) .* dt .>random));
#   v_temp[bond_events] .= true;
#   unbond_events = findall(x->x==1, (getfield.(v.u, :k_off).*dt .>random));
#   v_temp[unbond_events] .= false;

#   for i = 1:1:v.n                               # if I use setproperty! it doesn't work! Ask for help!
#     setfield!(v.u[i], :state, v_temp[i])
#   end

# end


