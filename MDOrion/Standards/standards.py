# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
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

from MDOrion.Standards.utils import ParmedData, MDStateData

from datarecord import OEPrimaryMolField

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField)


# ------------ Stage Standard Names ------------- #


class MDStageTypes:
    SETUP = 'SETUP'
    MINIMIZATION = 'MINIMIZATION'
    NVT = 'NVT'
    NPT = 'NPT'
    FEC = 'FEC'


# ------------ MD Engines ------------- #

class MDEngines:
    OpenMM = 'OpenMM'
    Gromacs = 'Gromacs'
    all = [OpenMM, Gromacs]


# ---------------- File  Name Standards -------------- #

class MDFileNames:
    topology = 'topology.oeb'
    state = 'state.pickle'
    trajectory = "trajectory.tar.gz"
    trajectory_conformers = "trajectory_confs.oeb"
    mddata = "data.tar.gz"

# ---------------- Field Standards -------------- #


class Fields:

    # The Title field is used to set the system name
    title = OEField("Title_OPLMD", Types.String, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The ID field should be used as identification number for ligands, proteins or complexes
    id = OEField("ID_OPLMD", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The SysID field is used to keep track of the system order
    sysid = OEField("SysID_OPLMD", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The ConfID field is used to identify a particular conformer
    confid = OEField("ConfID_OPLMD", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The Ligand field should be used to save in a record a ligand as an OEMolecule
    ligand = OEField("Ligand_OPLMD", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Ligand))

    # The ligand name
    ligand_name = OEField("Ligand_name_OPLMD", Types.String)

    # The Protein field should be used to save in a record a Protein as an OEMolecule
    protein = OEField("Protein_OPLMD", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Protein))

    # Primary Molecule
    primary_molecule = OEPrimaryMolField()

    # Parmed Structure, Trajectory, MDData and Protein trajectory conformers Fields
    if in_orion():
        pmd_structure = OEField('Structure_Parmed_OPLMD', Types.Int)
        trajectory = OEField("Trajectory_OPLMD", Types.Int)
        mddata = OEField("MDData_OPLMD", Types.Int)
        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Int)
    else:
        pmd_structure = OEField('Structure_Parmed_OPLMD', ParmedData)
        trajectory = OEField("Trajectory_OPLMD", Types.String)
        mddata = OEField("MDData_OPLMD", Types.String)
        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Chem.Mol)

    # The Stage Name
    stage_name = OEField('Stage_name_OPLMD', Types.String)

    # The Stage Type
    stage_type = OEField('Stage_type_OPLMD', Types.String)

    # Topology Field
    topology = OEField('Topology_OPLMD', Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.PrimaryMol))

    # Log Info
    log_data = OEField('Log_data_OPLMD', Types.String)

    # MD State
    md_state = OEField("MDState_OPLMD", MDStateData)

    # Collection is used to offload data from the record which mush be < 100Mb
    collection = OEField("Collection_ID_OPLMD", Types.Int)

    # Stage list Field
    md_stages = OEField("MDStages_OPLMD", Types.RecordVec)

    # Analysis Fields
    free_energy = OEField('FE_OPLMD', Types.Float,
                          meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

    free_energy_err = OEField('FE_Error_OPLMD', Types.Float,
                              meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

    floe_report = OEField('Floe_report_OPLMD', Types.String)

    floe_report_svg_lig_depiction = OEField("Floe_report_lig_svg_OPLMD", Types.String,
                                            meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))

    floe_report_label = OEField('Floe_report_label_OPLMD', Types.String)

