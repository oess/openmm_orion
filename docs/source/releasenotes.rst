#############
Release Notes
#############


v0.9.3
======================

General Notice
--------------------------------------------------------------------------------
* Upgrades to the trajectory analysis in floe ``Short Trajectory MD with Analysis`` expose key results to the Analyze page and 3D visualization page in Orion.
* Fixed bug in setting restraints in GROMACS cubes.
* Fixed bug in ordering results in the floe report in floe ``Short Trajectory MD with Analysis``.

New Floes
--------------------------------------------------------------------------------
* None


New Cubes
--------------------------------------------------------------------------------
* None


Cube Updates
--------------------------------------------------------------------------------
* Exposed MMPBSA ensemble average and standard deviation in the :ref:`cube_TrajPBSACube` so that it can be displayed in the Analyze page in Orion.
* :ref:`cube_MDTrajAnalysisClusterReport` now generates trajectory average and median molecules for protein and ligand, with one conformer for each major cluster. These are exposed int the 3D visualization page in Orion.
* :ref:`cube_ClusterOETrajCube` now exposes a link to the per-ligand floe report page so it is available in the Analyze page in Orion.
* In :ref:`cube_ComplexPrepCube` traditional references to the full periodic supermolecular ensemble as a "system" have been replaced with references to a "well" by analogy with an assay well.
* In :ref:`cube_MDFloeReportCube` the floe report now generates tiled links to individual ligands in the same order as the initial list of ligands.
* In :ref:`cube_MDFloeReportCube` the floe report tiles now show how many major clusters were found for each ligand.
* In :ref:`cube_MDNptCube` and :ref:`cube_MDNvtCube` the restraints are now correctly set in GROMACS for proteins consisting of multiple chains.


.. _2019.Oct: https://docs.eyesopen.com/toolkits/python/releasenotes/releasenotes/index.html
.. _OpenEye Toolkits: https://docs.eyesopen.com/toolkits/python/index.html