#############
Introduction
#############

This package is for running MD simulations useful for drug discovery,
with floes geared towards MD experts (General MD)
and floes also suitable for MD non-experts (Specialized MD).
Aside from the most general floes, the underlying design of the
floes is towards simulating complexes between a protein and a set
of small-molecule ligands in explicit water, with periodic boundary
conditions. The solutes are assumed to be prepared and ready for
such physics-based modeling.

Specialized MD
--------------
These floes carry out more focused MD tasks frequently needed in drug discovery,
particularly in lead optimization. 
Suitable for either MD experts or non-experts, the
protein and ligand inputs do 
require preparation appropriate for physics-based modeling but less
stringent than for General MD. These floes perform typical analyses of
the results suitable for quick, actionable interpretation by non-experts
in addition to access to the MD trajectory for more specialized
interpretation by experts.

General MD
--------------
These floes have the usual
stringent requirements for the preparation for input datasets
and produce simply the MD trajectory and standard output from the MD engine.
This places a relatively high onus on the user to prepare MD-ready input
and to access and analyze
the results using their own procedures.

Important notes about using this package
===================================
* Editing these floes in the Orion UI will probably not work due to the
  complex data structure used in the data stream between cubes.
