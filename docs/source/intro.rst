#############
Introduction
#############

This package is for running MD simulations useful for drug discovery,
enabling both MD experts and MD non-experts.
Floes for MD non-experts are in the category "Specialized MD": these
carry out more focused MD tasks frequently needed in drug discovery,
particularly in lead optimization. The protein and ligand inputs
require preparation ready for physics-based modeling but less
stringent than for General MD. These floes perform typical analyses of
the results suitable for quick, actionable interpretation by non-experts
in addition to access to the MD trajectory for more specialized
interpretation by experts.
Floes for MD experts are in the category "General MD": they have
stringent requirements for the preparation for input datasets
and produce simply the MD trajectory and standard output from the MD engine.
This places a relatively high onus on the user to prepare MD-ready input
and to access and analyze
the results using their own procedure.

Important notes about using this package
===================================
* Editing these floes in the Orion UI will probably not work due to the
  complex data structure used in the data stream between cubes.