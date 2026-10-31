---
title: 'CellAdhesion: A Julia package for modelling of stochastic interfaces'
tags:
  - Julia
  - adhesion
  - stochastic simulation
  - bond kinetics
  - mechanobiology
authors:
  - name: Alessandra Bonfanti
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Eleni Papafilippou
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Alexandre Kabla
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Engineering, University of Cambridge, United Kingdom
    index: 1
date: 31 October 2026
bibliography: paper.bib
---

# Summary

CellAdhesion is a Julia software package for modeling stochastic binding–unbinding dynamics across interacting interfaces. It provides a flexible computational framework for simulating how populations of bonds sharing a mechanical load evolve and rupture, under time-dependent loading. This makes the software broadly applicable to problems involving adhesive interactions, molecular linkages, or any system governed by probabilistic attachment–detachment processes.

The package includes multiple mechanistic models for the force dependent binding and unbinding rates, enabling users to evaluate alternative hypotheses and explore parameter regimes within a unified interface. It supports heterogeneous bond populations (allowing junctions composed of elements with distinct mechanical or kinetic properties) and bond clustering, where cooperative effects play a key role. These capabilities allow users to examine interface behavior from the single-bond scale to collective, emergent phenomena.

By combining performance, extensibility, and transparent model implementations, CellAdhesion offers a general-purpose platform for studying stochastic interface dynamics in biological, and engineered systems.

# Statement of need

Many scientific and engineering problems involve interfaces whose behaviour is governed by the stochastic formation and rupture of microscopic bonds. Modelling of such adhesive interactions play central roles in systems where probabilistic attachment-detachment events determine macroscopic behaviors, such as in biological adhesion [REF NEEDED] and self-healing polymers [REF NEEDED].

Despite their broad relevance, computational tools for simulating these stochastic dynamics remain fragmented. Early theoretical work, such as Bell's seminal study [REF NEEDED], focused on the binding and unbinding of a uniform population of molecular bonds under a constant load, emphasizing the stability and lifetime of such populations and allowing analytical solutions in idealized scenarios. More recent efforts have extended these models to capture complex behaviors, including heterogeneous bond populations, clustering effects, mixtures of bond types, and dynamic or non-uniform loading conditions, often relying on stochastic simulation techniques. However, implementations are frequently ad-hoc, embedded in domain-specific code, or tailored to a single class of problems, limiting their general applicability. A single, integrated computational framework that combines these capabilities would provide a versatile platform for systematically addressing a wide range of open questions in adhesive interactions.

Existing stochastic simulation libraries provide general Gillespie-type methods, but lack domain-specific support for force-dependent unbinding, heterogeneous bond populations, or mechanically coupled junction evolution. Conversely, particle-based or molecular-dynamics tools can represent bond breakage, yet are computationally too expensive for simulating large adhesion clusters or for conducting ensemble-level lifetime studies. As a result, researchers often resort to ad-hoc scripts or highly specialized in-house code, limiting reproducibility, extensibility, and cross-domain applicability.

A dedicated framework (Figure 1) is therefore needed in order to:

- construct adhesion interfaces with diverse spatial organizations, such as uniform, or clustered distributions;
- assign heterogeneous kinetic properties to different bonds, capturing the simultaneous existence of different molecular species, capturing distinct binding and unbinding dynamics within the same junction; and
- apply any time-dependent loading profile (e.g. constant, ramp).

Such capabilities are essential for computing the full dynamical evolution of an adhesive junction at the microscopic scale, including its ultimate lifetime under mechanical load.

![High level schematic of the workflow of a stochastic interface simulation.](figures/JOSS_CellAdhesion.png){ width=50% }

# Implementation

The library is implemented in the Julia programming language and follows a modular design that cleanly separates data structures, physical models, and numerical simulation routines. All functionality is organized within the `CellAdhesion` module.

**Data structures**

The fundamental unit of the implementation is the `Bond` type, which stores the mechanical state of an individual adhesion molecule (open/closed and applied force) together with its associated kinetic model. To represent systems containing multiple adhesion types, the library defines an abstract `BondModel` type from which specific kinetic models (e.g. slip bonds) inherit. This structure allows heterogeneous bonds to be grouped while preserving type-specific parameters.

Interfaces are represented using the `Cluster` type, which holds a dynamic collection of `Bond` instances and tracks the corresponding cluster-level force and state. Larger interfaces can be composed of multiple `Cluster` objects. A schematic summary of these data structures and their relationships is shown in Figure 2.

**Loading force**

The implementation supports arbitrary external loading conditions, including constant loads and fully user-defined time-dependent force profiles. At each simulation step, the current load is evaluated and passed to the corresponding bond models to compute instantaneous force-dependent off-rates. Within a cluster, the applied force can be redistributed among bonds either uniformly or through state-dependent rules that account for the configuration of neighboring open and closed bonds.

**Bond model kinetics**

Force-dependent binding and unbinding dynamics are implemented via model-specific rate functions. The library currently includes Bell-type slip-bond formulations but is designed to be readily extensible: users can introduce new models by defining additional `BondModel` subtypes with custom association and disassociation rate laws tailored to their physical system.

**Simulation of interface lifetime**

The evolution of the adhesive interface is simulated through stochastic updates of the bond states, governed by Monte Carlo sampling of the force-dependent binding and unbinding probabilities. The simulation advances in time until the entire bond cluster dissociates or a user-defined terminal time is reached, thereby providing both the temporal evolution of the junction and the resulting interface lifetime together with the corresponding rupture force.

![Schematic of the data structures implemented within the library.](figures/JOSS_data_structure.png){ width=75% }

# Research use

CellAdhesion has been used in recent research works, including fracture of biological tissues [@duque2024rupture] and resilience of biological tissues under cyclic loading [@papafilippoui2025]. Furthermore, it is in use for the study of cell-cell junction interfaces to assess the role of molecules heterogeneity in junction strength.

# Acknowledgements

A.B. was supported by the European Research Executive Agency via the MSCA Postdoctoral Fellowships (BIOFRAC, Project number: 101105376). E.P. was supported by the Engineering and Physical Sciences Research Council Centre for Doctoral Training in Sensor Technologies for a Healthy and Sustainable Future [EP/S023046/1]. We thank the past and present members of the Kabla laboratory for insightful discussion.

# References
