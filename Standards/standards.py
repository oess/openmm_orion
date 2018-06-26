# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from floe.api.orion import in_orion

from Standards.utils import ParmedData

from datarecord import OEPrimaryMolField

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

# from cuberecord.record_types.cube_types import OrionTypes
from cuberecord import Types as TypesCR


# ------------ Stage Standard Names ------------- #

class MDStageNames:
    SETUP = 'SETUP'
    MINIMIZATION = 'MINIMIZATION'
    NVT = 'NVT'
    NPT = 'NPT'
    FEC = 'FEC'

# ---------------- Field Standards -------------- #


class Fields:
    # The Title field is used to set the system name
    title = OEField("Title_OPLMD", Types.String, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The ID field should be used as identification number for ligands, proteins or complexes
    id = OEField("ID_OPLMD", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The Ligand field should be used to save in a record a ligand as an OEMolecule
    ligand = OEField("Ligand_OPLMD", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Ligand))

    # The Protein field should be used to save in a record a Protein as an OEMolecule
    protein = OEField("Protein_OPLMD", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Protein))

    # Primary Molecule
    primary_molecule = OEPrimaryMolField()

    # Parmed Structure Field
    structure = OEField('Structure_Parmed_OPLMD', ParmedData)

    # The Stage Name
    stage_name = OEField('Stage_name_OPLMD', Types.String)

    # Topology Field
    topology = OEField('Topology_OPLMD', Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.PrimaryMol))

    # Log Info
    log_data = OEField('Log_data_OPLMD', Types.String)

    # MD System Field
    md_system = OEField("MDSystem_OPLMD", Types.Record)

    md_system.get_type()

    # Trajectory
    if in_orion():
        trajectory = OEField("Trajectory_OPLMD", TypesCR.Orion.LargeFile)
    else:
        trajectory = OEField("Trajectory_OPLMD", Types.String)

    # Stage list Field
    md_stages = OEField("MDStages_OPLMD", Types.RecordVec)

    # Stage Field
    md_stage = OEField("MDStages_OPLMD", Types.Record)

    yank_analysis = OEField("Yank_Analysis_OPLMD", Types.String)


# ---------------- Record Standards -------------- #

# The MDSystemRecord class holds the system topology as an OEMol and the system
# parametrization by using a Parmed Structure object

class MDRecords:
    class MDSystemRecord(OERecord):

        def __init__(self, molecule, structure):
            super().__init__()
            self.set_value(Fields.topology, molecule)
            self.set_value(Fields.structure, structure)

    class MDStageRecord(OERecord):

        def __init__(self, name, system_record, log=None, trajectory=None):
            super().__init__()
            self.set_value(Fields.stage_name, name)
            self.set_value(Fields.md_system, system_record)
            if log is not None:
                self.set_value(Fields.log_data, log)
            if trajectory is not None:
                self.set_value(Fields.trajectory, trajectory)
