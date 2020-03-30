The Plain MD protocol performs MD simulations given one or more
complete molecular systems as input, each to be treated in its entirety as a solute.
The solute need to have coordinates, all atoms, and correct chemistry.
Each molecular system can have multiple conformers but each conformer will be
run separately as a different solute.
Proteins need to be prepared to an MD standard: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
The input system is solvated and parametrized according to the
selected force fields. A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up and equilibration stages positional harmonic restraints are
applied. At the end of the equilibration stages a short
(default 2ns) production run is performed on the unrestrained system.
