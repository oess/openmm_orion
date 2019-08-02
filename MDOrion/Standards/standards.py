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


# Orion Hidden meta data options
_metaHidden = OEFieldMeta(options=[Meta.Display.Hidden])
_metaIDHidden = OEFieldMeta(options=[Meta.Source.ID, Meta.Display.Hidden])
_metaProtHidden = OEFieldMeta(options=[Meta.Hints.Chem.Protein, Meta.Display.Hidden])


# ---------------- Field Standards -------------- #
class Fields:

    # The Title field is a string name for the well which used to compose file names
    title = OEField("Title_OPLMD", Types.String, meta=_metaIDHidden)

    # The wellid field is a unique integer for each well (final system for simulation)
    wellid = OEField("WellID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The ligid field is a unique integer used to keep track of the ligand input order
    ligid = OEField("LigID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The ConfID field is used to identify a particular conformer
    confid = OEField("ConfID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The Ligand field should be used to save in a record a ligand as an OEMolecule
    ligand = OEField("Ligand_OPLMD", Types.Chem.Mol, meta=OEFieldMeta(options=[Meta.Hints.Chem.Ligand,
                                                                               Meta.Display.Hidden]))

    # The ligand name
    ligand_name = OEField("Ligand_name_OPLMD", Types.String, meta=_metaHidden)

    # The protein field should be used to save in a record a Protein as an OEMolecule
    protein = OEField("Protein_OPLMD", Types.Chem.Mol, meta=_metaProtHidden)

    # The protein name
    protein_name = OEField("Protein_name_OPLMD", Types.String, meta=_metaHidden)

    # The super-molecule for the entire Well (ie the final system for simulation)
    well = OEField("Well_OPLMD", Types.Chem.Mol, meta=_metaHidden)

    # Primary Molecule
    primary_molecule = OEPrimaryMolField()

    # Parmed Structure, Trajectory, MDData and Protein trajectory conformers Fields
    if in_orion():
        pmd_structure = OEField('Structure_Parmed_OPLMD', Types.Int, meta=_metaHidden)
        trajectory = OEField("Trajectory_OPLMD", Types.Int, meta=_metaHidden)
        mddata = OEField("MDData_OPLMD", Types.Int, meta=_metaHidden)
        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Int, meta=_metaHidden)
    else:
        pmd_structure = OEField('Structure_Parmed_OPLMD', ParmedData, meta=_metaHidden)
        trajectory = OEField("Trajectory_OPLMD", Types.String, meta=_metaHidden)
        mddata = OEField("MDData_OPLMD", Types.String, meta=_metaHidden)
        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Chem.Mol, meta=_metaHidden)

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

    # Collection is used to offload data from the record which must be < 100Mb
    collection = OEField("Collection_ID_OPLMD", Types.Int, meta=_metaHidden)

    # Stage list Field
    md_stages = OEField("MDStages_OPLMD", Types.RecordVec, meta=_metaHidden)

    # Analysis Fields
    free_energy = OEField('FE_OPLMD', Types.Float,
                          meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

    free_energy_err = OEField('FE_Error_OPLMD', Types.Float,
                              meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

    floe_report = OEField('Floe_report_OPLMD', Types.String, meta=_metaHidden)

    floe_report_svg_lig_depiction = OEField("Floe_report_lig_svg_OPLMD", Types.String,
                                            meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))

    floe_report_label = OEField('Floe_report_label_OPLMD', Types.String, meta=_metaHidden)

    floe_report_URL = OEField('Floe_report_URL_OPLMD', Types.String, meta=OEFieldMeta(options=[Meta.Hints.URL]))

    class Analysis:

        # The OETraj Field is for the record containing Traj OEMols and energies
        oetraj_rec = OEField("OETraj", Types.Record, meta=_metaHidden)

        # The TrajIntE Field is for the record containing Traj interaction energies
        oeintE_rec = OEField("TrajIntE", Types.Record, meta=_metaHidden)

        # The TrajPBSA Field is for the record containing Traj PBSA energies
        oepbsa_rec = OEField("TrajPBSA", Types.Record, meta=_metaHidden)

        # The TrajClus Field is for the record containing Traj ligand clustering results
        oeclus_rec = OEField("TrajClus", Types.Record, meta=_metaHidden)

        # The AnalysesDone Field is for a list of the analyses that have been done
        analysesDone = OEField("AnalysesDone", Types.StringVec, meta=_metaHidden)

        # The Lig_Conf_Data Field is for the record containing Traj conf data for all confs
        oetrajconf_rec = OEField("Lig_Conf_Data", Types.RecordVec, meta=_metaHidden)

        # The vector of ligand Traj RMSDs from the initial pose
        lig_traj_rmsd = OEField('LigTrajRMSD', Types.FloatVec,
                                   meta=OEFieldMeta().set_option(Meta.Units.Length.Ang))

        # mmpbsa ensemble average over the whole trajectory
        mmpbsa_traj_mean = OEField('MMPBSATrajMean', Types.Float,
                                   meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

        metaMMPBSA_traj_std = OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol)
        metaMMPBSA_traj_std.add_relation(Meta.Relations.ErrorsFor, mmpbsa_traj_mean)
        mmpbsa_traj_std = OEField('MMPBSATrajStdev', Types.Float, meta=metaMMPBSA_traj_std)

        # The number of major clusters found
        n_major_clusters = OEField("n major clusters", Types.Int)

        # Trajectory cluster averages and medians of protein and ligand
        ClusLigAvg_fld = OEField('ClusLigAvgMol', Types.Chem.Mol)
        ClusProtAvg_fld = OEField('ClusProtAvgMol', Types.Chem.Mol)
        ClusLigMed_fld = OEField('ClusLigMedMol', Types.Chem.Mol)
        ClusProtMed_fld = OEField('ClusProtMedMol', Types.Chem.Mol)


def get_meta_attributes(record, field_name):
    field_with_meta = record.get_field(field_name, include_meta=True)
    meta_from_field = field_with_meta.get_meta()
    meta_dict = meta_from_field.to_dict()
    return meta_dict
