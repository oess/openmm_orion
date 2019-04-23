#############
Release Notes
#############


v0.8.3
======================

General Notice
--------------------------------------------------------------------------------
* Upgrades to ``OpenEye-orionplatform==0.1.14``

New Floes
--------------------------------------------------------------------------------
* Upgrades to ``OpenEye-orionplatform==0.1.14``


New Cubes
--------------------------------------------------------------------------------
* New :ref:`cube_ConformerSplitterCube` and :ref:`cube_ConformerMergerCube` cubes which preserve fields added during
  calculations on the split conformers.
    * previous cubes for Splitting and Merging conformers which did not preserve fields have been removed.


Cube Updates
--------------------------------------------------------------------------------
* Added new parameters `mode` and `rotor_offset` in the :ref:`cube_OmegaCube`.
* Fixed a parameter mismatch issue with the :ref:`cube_AssignChargesCube`.
* The :ref:`cube_TorsionScanCube` is updated to use a single tag to identify a torsion, instead of four tags. It also
now has options to generate conformations using Omega prior to torsion driving or use a provided conformer ensemble. By
default the cube reports both the MMFF/Sheffield energy as well as the MMFF/PB energy for each output conformer.
* Fixed :ref:`cube_ZapLigandPrepCube` to fail gracefully when it cannot prep a molecule.
* The :ref:`cube_TorsionScanCube` is updated to use a signle tag to identify a torsion, instead of four tags.
* New parameters `in_strain_mol_field` and `out_strain_mol_field` added to the :ref:`cube_FreeformStrainCube` and the
:ref:`cube_FreeformEneStrainCube` cubes to enable using separate input conformers for freeform free energy and strain calculations.
* New parameter `use_inp_ens` have been added, to enable using user provided conformer ensemble in calculations, to the
following freeform cubes:

  * :ref:`cube_FreeformConfEntropyCube`
  * :ref:`cube_FreeformConfFreeEnergiesCube`
  * :ref:`cube_FreeformStrainCube`
  * :ref:`cube_FreeformEneStrainCube`
  * :ref:`cube_FreeFormConfPrepEnsembleCube`
  * :ref:`cube_FreeFormConfPreOptimizeCube`
  * :ref:`cube_FreeFormConfRemoveDupsCube`
  * :ref:`cube_FreeFormConfOptimizeCube`
  * :ref:`cube_FreeFormConfEnergiesCube`


.. _2018.Oct: https://docs.eyesopen.com/toolkits/python/releasenotes/releasenotes/index.html
.. _OpenEye Toolkits: https://docs.eyesopen.com/toolkits/python/index.html