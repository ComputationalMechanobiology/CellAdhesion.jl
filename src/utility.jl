
export print_cluster, plot_cluster, plot_force_distribution


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





"""
  plot_force_distribution(v, p, x, y, cmap)

  Plot bond tensions for a cluster, using the colormap `cmap` for closed bonds
  and gray squares for open bonds.

"""
# Leaf cluster: bonds stacked along y, colored by tension.
function plot_force_distribution!(v::Cluster{Bond{T}}, p, x, y = 0.0; cmap = :viridis,
                       show_colorbar::Bool = true) where T <: BondModel
  states, forces = bond_state_force(v)
  y_coords = y .+ (0:(v.n - 1)) .* v.l
  closed     = findall(states)
  open_bonds = findall(.!states)

  unit_square(x, y) = Shape(
      [x - 0.5, x + 0.5, x + 0.5, x - 0.5],
      [y - 0.5, y - 0.5, y + 0.5, y + 0.5]
  )

  if !isempty(closed)
      for i in closed
          plot!(p, unit_square(x, y_coords[i]);
                fill_z = forces[i],
                seriescolor = cmap,
                linecolor = :transparent,
                colorbar = show_colorbar,
                colorbar_title = "Tension",
                label = "")
      end
  end

  if !isempty(open_bonds)
      for i in open_bonds
          plot!(p, unit_square(x, y_coords[i]);
                fillcolor = "#e0e0e0",
                linecolor = :transparent,
                label = "")
      end
  end
end

# Hierarchical cluster: sub-clusters laid out horizontally (different x).
function plot_force_distribution!(v::Cluster, p, x, y ; cmap = :viridis, show_colorbar::Bool = true)
  running_y = 0     
  for i = 1:v.n
    plot_force_distribution!(v.u[i], p, x , y + running_y; cmap = cmap, show_colorbar = show_colorbar && i == 1)
    running_y += v.u[i].n * v.u[i].l + v.l
  end
end

# Wrapper that creates the plot and returns it.
function plot_force_distribution(v, p, x, y, cmap)
  p = plot!(p, legend = false, framestyle = :box)
  plot_force_distribution!(v, p, x, y; cmap = cmap)
  xlabel!(p, "time")
  ylabel!(p, "Bond index")
  return p
end