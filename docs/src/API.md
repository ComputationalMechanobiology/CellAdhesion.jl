# API
```@meta
EditURL = "https://github.com/ComputationalMechanobiology/CellAdhesion.jl/edit/main/docs/src/API.md"
```

## Main CellAdhesion data structures
```@docs
Bond
Cluster
```

## Processing functions
```@docs
CellAdhesion.update!
CellAdhesion.runcluster
CellAdhesion.bond_state_force
CellAdhesion.save_cluster_to_json
CellAdhesion.load_from_json
```

## Utility functions
```@docs
CellAdhesion.state!
CellAdhesion.print_cluster
```

## Serialization internals
```@docs
CellAdhesion._save_bond_to_dict
CellAdhesion._save_cluster_to_dict
CellAdhesion._load_bond_from_dict
CellAdhesion._load_cluster_from_dict
```

## Dynamics functions
```@docs
CellAdhesion.setforce!
CellAdhesion.force_global
CellAdhesion.force_local
CellAdhesion.distance
```




