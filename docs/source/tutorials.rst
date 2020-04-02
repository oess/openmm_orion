#############
Tutorials
#############

These are tutorials for how to use specific floes.

Short Trajectory MD with Analysis
=================================

The Floe *Short Trajectory MD with Analysis* (STMD) is for pose validation.
You should start with a ligand already well posed in the active site,
with correct chemistry (including formal charges)
and all atoms (including explicit hydrogens, which determines a
specific tautomer). As such, this pose represents an hypothesis for
ligand binding, with protein-ligand interactions in place, which will
be evaluated by this floe to be validated (or not) as a good pose.
This floe is not primarily asking the question "what is the pose?"
but rather "is this a good pose?", seeking an answer by running a
short MD trajectory and assessing the ligand pose especially by
comparison to the starting (input) pose.

Ligand Input
------------

Just to be able to run, this floe requires ligands to have
reasonable 3D coordinates, all atoms, and correct chemistry
(in particular bond orders and formal charges).
If the ligands already have good atomic partial charges
(we recommend RESP or AM1-BCC_ELF10 charges),
we recommend using these for STMD as opposed to re-charging
them in the STMD floe.
Given that this floe only runs a very short timescale (default 2 ns)
to evaluate the input pose,
it is preferable that the input pose be well refined.
Although bad clashes
(or poor positioning for interactions which you know need to be there)
can be (and often are) cleaned up by even this short trajectory,
it starts off the "evaluation" purpose of the floe on the wrong foot
by giving a poor comparator.
Poor initial poses might even be considered outside the scope of this floe
given how short is the default timescale.
This is why we **strongly recommend** that docked poses be
subsequently minimized in the active site **before input to STMD**.
This will resolve high gradients
(usually clashes) with the protein and to allow protein-ligand
interactions to optimize in the context of a good force field.
It is possible that even with this, docked-pose starting points
could be triaged prior to the extra effort and expense of STMD.

Protein Input
-------------
All the MD floes require correctly prepared protein up to "MD ready" standards.
This begins with the normal prerequisites for physics-based modeling:
protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved or capped.
Of course, protein side chain formal charges and protonation
at this point determine their tautomeric state.
Additionally, cofactors and structured internal waters are also important to include,
not only those in the immediate vicinity of the ligand and active site
but also distally because they can have an important effect on the
protein structure and dynamics over the course of the MD.
We **strongly recommend** using an appropriate *Spruce* floe for protein preparation.

.. warning::

   Unfortunately, proteins with covalently bound ligands or cofactors are currently not tractable

How to use this floe
--------------------
The structure of the STMD floe is shown in figure
`Structure of the STMD floe`.
Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer separately,
and the complex is solvated and parametrized according to
the selected force fields.
A minimization stage is performed on the system followed by
a warm up (NVT ensemble) and three equilibration stages (NPT ensemble).
In the minimization, warm up, and equilibration stages,
positional harmonic restraints are applied on the ligand and protein.
At the end of the equilibration stages a short (default 2ns) production run
is performed on the unrestrained system.
The production run is then analyzed in terms of interactions between
the ligand and the active site and in terms of ligand RMSD,
after fitting the trajectory based on active site C_alphas.

.. figure_STMD_floe:

.. figure:: ./images/STMD_floe.png
   :width: 800px
   :align: center
   :alt: Structure of the STMD floe

   **Structure of the STMD floe**

After selecting the *Short Trajectory MD with Analysis* floe in the Orion UI,
you will be presented with a job form with parameters to select.
Aside from the essential user-defined parameters relating to jobname,
input (protein and ligand datasets as described above), and
output (output and failure dataset names),
all other parameters have defaults preset to values that
we feel are suitable for the non-expert MD user
(or expert user for that matter) so launching the floe at this point is reasonable.
Of the other top-level parameters a few will be discussed here:

    * Protein_title (no default): Here is where you can put a handy short name for the protein to use in molecule titles (e.g. "Bace" instead of "beta-secretase").

    * Charge_ligands (default *True*): If your input ligands already have good atomic partial charges (e.g. `RESP` or `AM1-BCC_ELF10`), set this to *False* to have the floe use the existing ligand charges.

    * Ligand_forcefield (default *OpenFF1.0.0*): This forcefield choice has a strong impact on the results. We recommend the most recent version of the OpenFF force field from the *Open Force Field Initiative*.

    * Md_engine (default *OpenMM*): Gromacs is the other alternative but we recommend OpenMM because HMR works with it but not with Gromacs.

    * Hmr: Hydrogen Mass Repartitioning (HMR) gives a two-fold speedup and reduces cost. We recommend leaving it on.

We make the other top-level parameters available for expert users.
