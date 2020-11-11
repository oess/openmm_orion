The Short Trajectory MD (STMD) protocol performs MD simulations given
a prepared protein and a set of posed and prepared ligands as input.
The ligands need to have reasonable 3D coordinates, all atoms, and correct chemistry
(in particular bond orders and formal charges). Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand. The starting poses should not have very high gradients,
in particular no bad clashes with the protein.
The protein needs to be prepared to MD standards: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved or capped. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the complex is solvated and parametrized according to the
selected force fields. A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up, and equilibration stages, positional harmonic restraints are
applied on the ligand and protein. At the end of the equilibration stages a short
(default 2ns) production run is performed on the unrestrained system.
The production run is then analyzed in terms of interactions between the
ligand and the active site and in terms of ligand RMSD after fitting the trajectory
based on active site C_alphas.
