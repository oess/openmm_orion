The Protein-Ligand MD floe performs MD simulations given
a prepared protein and a set of posed and prepared ligands as input.
The ligands need to have coordinates, all atoms, and correct chemistry. Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand.
The protein needs to be prepared to MD standards: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the complex is solvated and parametrized according to the
selected force fields. A minimization stage is peformed on the system followed
by a warm up (NVT ensemble) and several equilibration stages (NPT ensemble). In the
minimization, warm up, and equilibration stages, positional harmonic restraints are
applied on the ligand and protein. At the end of the equilibration stages a 
production run (by default only 2 ns) is performed on the unrestrained system.
The trajectory and final state are written out in the results record;
no analysis is performed.

Required Input Parameters:
* A ligand dataset of prepared ligands posed in the protein active site.
* A protein dataset of the prepared MD-ready protein structure, including cofactors and structured waters.
