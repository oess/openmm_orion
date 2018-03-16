# Orion OpenMM data records v0.1

---------------------------

To ensure interoperability of OpenMM-based cubes and floes making use of `OERecord` to exchange information, we propose a common set of `OERecord` and `OEField` definitions that can be imported from a common library via
```python
from openmm_orion.standards import *
```
When possible, we encourage all users building OpenMM-based floes to utilize these records and fields, and to avoid naming collisions to ensure interoperability.

## Fields (`OEField`)

The following `OEField` objects are defined for the storage of simulation data:
* `topology_field = OEField('topology', Types.OEMol)` : an `OEMol` object describing the chemical matter contained in this system, and may contain multiple molecules (such as protein and solvent)
* `structure_field = OEField('structure', parmed.Structure)`: a ParmEd structure describing the positions of all atoms, box vectors, and forcefield parameters for a corresponding topology
* `logdata_field = OEField('logdata', Types.TextFile)` : text log data
* `stage_name_field = OEField('stage_name', Types.Text)` : a text field describing the current stage of preparation or simulation

**QUESTION:** Do we want to standardize the potential naming scheme for stages?

## Records (`OERecord`)

### `MDSystemRecord`

An `MDSystemRecord` contains the information necessary to run (or extend) a molecular dynamics simulation.
The record consists of the following fields:
* `topology_field` (field name: `topology`): The chemical matter contained in the system, which may contain multiple molecules
* `structure_field` (field name: `structure`): The positions, box vectors, and forcefield parameters for the system

The `structure_field` and `topology_field` must have the same number of atoms, atom connectivity, and atomic positions.
If virtual sites are in use, the `parmed.Structure` object of the `structure_field` will have more particles than the number of atoms in the `topology_field`, as the `OEMol` will represent only the chemical atoms in the system.

**QUESTION:** Should we allow the `parmed.Structure` object to lack parameters and still be a valid `MDSystemRecord`?
This would allow parameterization cubes to both accept and return an `MDSystemRecord`, but it would mean that we aren't guaranteed to be able to simulate the resulting `MDSystemRecord`.

### `MDStageRecord`

An `MDStageRecord` is the result of a simulation that produces a result and is capable of being extended.
The record consists of the following fields:
* `stage_name_field` (field name: `stage_name`): The name of the stage that produced this result
* (field_name: `trajectory`) : An `OELargeFile` containing the trajectory output data in an [MDTraj](http://mdtraj.org/)-readable format. We recommend use of the [MDTraj HDF5 format](http://mdtraj.org/1.9.0/hdf5_format.html) since it provides the ability to [write multiple properties to the HDF5 file](http://mdtraj.org/1.9.0/api/generated/mdtraj.formats.HDF5TrajectoryFile.html#mdtraj.formats.HDF5TrajectoryFile.write).
* `log_data_field` (field name: `log_data`): The text logfile produced by the stage
* `MDSystemRecord` (field name: `system_record`): The `MDSystemRecord` representing the final state of the computation

## Notes

The following diagram illustrates these concepts:

[![Data Record MD](https://github.com/oess/openmm_orion/tree/gcalabro_data_record/MD_Data_Record/images/Plan_MD_DataRecord.png)](https://github.com/oess/openmm_orion/tree/gcalabro_data_record/MD_Data_Record/images/Plan_MD_DataRecord.png)
