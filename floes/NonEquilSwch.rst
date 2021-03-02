The *Non-Equilibrium Swiching* Floe performs Relative Binding Free Energy (RBFE)
Calculations using the Non-Equilibrium Switching (NES) refined by the de Groot lab
(Gapsys et al., Chem. Sci., 2020, 11, 1140-1152). This floe takes three inputs:
(1) An Orion dataset containing an equilibrium run for each bound ligand participating in the NES run.
(2) An Orion dataset containing an equilibrium run for each unbound ligand participating in the NES run.
(3) A text file of the desired transformations of one ligand into another ("edges").
The file has one line per transformation, of format "ligA >> ligB"
where "ligA" and "ligB" are the ligand names for the ligands to be transformed. 
These names must correspond exactly to those in the unbound and bound ligand equilibration datasets.
The floe will draw a number of starting snapshots from the bound and unbound trajectories,
generate an RBFE alchemical transformation from ligA into ligB,
and carry out the NES fast transformation of ligA into ligB, and vice versa, for each of the snapshots.
The speed of the NES transformation and the number of snapshots transformed
can be adjusted from default values by the user at runtime.
The floe outputs a dataset containing a results record for each edge with the calculated
RBFE free energy difference (delta delta G) between the ligands for each edge.
The floe also outputs an HTML floe report having a top-level summary for each edge
linked to a more detailed description.
