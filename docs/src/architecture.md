# Architecture

```@meta
EditURL = "https://github.com/ComputationalMechanobiology/CellAdhesion.jl/edit/master/docs/src/architecture.md"
```

CellAdhesion is built around three main data types:
+ [`Bond`](@ref): it is a mutable structure that contains information regarding the state of the bond (open = false or closed = true), the force applied to the bond (f), the model that discribes the dynamics of the bond (binding-unbinding probabilities).
+ [`Cluster`](@ref): it is a mutable structure that contains a unit element vector (this can be of type Bond if the junction is made of one bond cluster, or Cluster if the junction is made of clusters), the state of the cluster (open = false, closed = true), the force applied to the Cluster (f), the model by which the force is distributed within bonds (f_model), the number of unit elements (n), the distance between unit elements (l). 

