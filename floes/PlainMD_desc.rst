The *Solvate and Run MD* Floe performs MD simulations given one or more
complete solutes as input, each to be treated in its entirety in a separate simulation.
Each solute is solvated in a
The solute need to have coordinates, all atoms, and correct chemistry.
Each solute can have multiple conformers but each conformer will again be
run separately as a different solute.
Solutes need to be prepared to an MD standard: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved or capped. Crystallographic internal waters should be retained where
possible. Non-protein components must likewise have all atoms, and correct chemistry
(in particular bond orders and formal charges).
The parametrization of some common nonstandard residues is partially supported.
The input system is solvated and parametrized according to the
selected force fields. A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up and equilibration stages positional harmonic restraints are
applied. At the end of the equilibration stages a
production run (default 2ns) is performed on the unrestrained system.
**The output produced by the floe is not compatible with the Analyze Protein-Ligand MD floe**
