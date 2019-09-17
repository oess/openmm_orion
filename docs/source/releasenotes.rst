#############
Release Notes
#############


v0.9.4
======================

General Notice
--------------------------------------------------------------------------------
* The MD Spruce Prep Floe has been removed. Proteins must be prepared with the Spruce Prep floes available in the Classic
   floes now.

* The calculation of MMPBSA can now be also performed by using explicit waters (still experimental)

New Floes
--------------------------------------------------------------------------------

* None

--------------------------------------------------------------------------------

Floe Updates
--------------------------------------------------------------------------------

* The MD Spruce Pre Floe has been removed

--------------------------------------------------------------------------------

New Cubes
--------------------------------------------------------------------------------
* A new cube has been developed to check the record size before to be written to the Orion backend
    avoiding floe failures. The new cube has been added to all the floes for sanity check.

* A new cube to estimate the water number around a ligand-protein complex has been developed. The cube is
    used in the MMPBSA calculation with the explicit water flag set on

Cube Updates
--------------------------------------------------------------------------------
* Hint interactions and Styles have been removed from receptors and ligands in the Protein, Ligand and FF parametrization
    setting cubes that could cause problems along the MD analysis stages (debugging is in progress)

* A bug has been fixed in the ligand Elf10 charging cube that was causing problems when carboxylic acid was present
    in a ligand to be charged

* The Trajectory to OEMol, Interaction Energies and PBSA calculation cubes have been updated to account for the explicit
    water in the new MMPBSA calculation

