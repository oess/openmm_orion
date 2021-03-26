# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
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

from MDOrion.Standards.utils import ParmedData, MDStateData, DesignUnit, MDComponentData

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


class MDStageNames:
    ForceField = "System Parametrization"
    Minimization = "System Minimization"
    WarmUp = "WarmUp"
    EquilibrationI = "EquilibrationI"
    EquilibrationII = "EquilibrationII"
    EquilibrationIII = "EquilibrationIII"
    Production = "Production"


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

# ---------------- Collection  Name Standards -------------- #


class CollectionsNames:
    none = ''
    md = 'MD_OPLMD'
    nes = 'NES_OPLMD'


# Orion Hidden meta data options
_metaHidden = OEFieldMeta(options=[Meta.Display.Hidden])
_metaIDHidden = OEFieldMeta(options=[Meta.Source.ID, Meta.Display.Hidden])
_metaProtHidden = OEFieldMeta(options=[Meta.Hints.Chem.Protein, Meta.Display.Hidden])


# ---------------- Field Standards -------------- #
class Fields:

    # The LigInitialRecord Field is for the initial ligand record read in at the start
    ligInit_rec = OEField("LigInitial", Types.Record, meta=_metaHidden)

    # The Title field is a string name for the flask which used to compose file names
    title = OEField("Title_OPLMD", Types.String, meta=_metaIDHidden)

    # The flaskid field is a unique integer for each flask (final system for simulation)
    flaskid = OEField("FlaskID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The ligid field is a unique integer used to keep track of the ligand input order
    ligid = OEField("LigID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The ConfID field is used to identify a particular conformer
    confid = OEField("ConfID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The Ligand field should be used to save in a record a ligand as an OEMolecule
    ligand = OEField("Ligand_OPLMD", Types.Chem.Mol, meta=OEFieldMeta(options=[Meta.Hints.Chem.Ligand, Meta.Display.Hidden]))

    # The ligand name
    ligand_name = OEField("Ligand_name_OPLMD", Types.String, meta=_metaHidden)

    # The protein field should be used to save in a record a Protein as an OEMolecule
    protein = OEField("Protein_OPLMD", Types.Chem.Mol, meta=_metaProtHidden)

    # The protein name
    protein_name = OEField("Protein_name_OPLMD", Types.String, meta=_metaHidden)

    # The super-molecule for the entire flask (ie the final system for simulation)
    flask = OEField("Flask_OPLMD", Types.Chem.Mol, meta=_metaHidden)

    # Primary Molecule
    primary_molecule = OEPrimaryMolField()

    # Parmed Structure, Trajectory, MDData and Protein trajectory conformers Fields
    if in_orion():
        pmd_structure = OEField('Structure_Parmed_OPLMD', Types.Int, meta=_metaHidden)
        trajectory = OEField("Trajectory_OPLMD", Types.Int, meta=_metaHidden)
        mddata = OEField("MDData_OPLMD", Types.Int, meta=_metaHidden)
        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Int, meta=_metaHidden)
        extra_data_tar = OEField("ExtraData_OPLMD", Types.Int, meta=_metaHidden)
    else:
        pmd_structure = OEField('Structure_Parmed_OPLMD', ParmedData, meta=_metaHidden)
        trajectory = OEField("Trajectory_OPLMD", Types.String, meta=_metaHidden)
        mddata = OEField("MDData_OPLMD", Types.String, meta=_metaHidden)
        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Chem.Mol, meta=_metaHidden)
        extra_data_tar = OEField("ExtraData_OPLMD", Types.String, meta=_metaHidden)

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

    # Design Unit Field
    design_unit = OEField('Design_Unit_OPLMD', DesignUnit)

    # Design Unit Field from Spruce
    # design_unit_from_spruce = OEField('du_single', Types.Blob)
    design_unit_from_spruce = OEField('designunit', Types.Chem.DesignUnit)

    # MD Components
    md_components = OEField('MDComponents_OPLMD', MDComponentData)

    # Collection is used to offload data from the record which must be < 100Mb
    # collection = OEField("Collection_ID_OPLMD", Types.Int, meta=_metaHidden)

    collections = OEField("Collections_ID_OPLMD", Types.JSONObject, meta=_metaHidden)

    # Stage list Field
    md_stages = OEField("MDStages_OPLMD", Types.RecordVec, meta=_metaHidden)

    floe_report = OEField('Floe_report_OPLMD', Types.String, meta=_metaHidden)

    floe_report_svg_lig_depiction = OEField("Floe_report_lig_svg_OPLMD", Types.String,
                                            meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))

    floe_report_label = OEField('Floe_report_label_OPLMD', Types.String, meta=_metaHidden)

    floe_report_URL = OEField('Floe_report_URL_OPLMD', Types.String, meta=OEFieldMeta(options=[Meta.Hints.URL]))

    floe_report_collection_id = OEField('Floe_report_ID_OPLMD', Types.Int, meta=_metaHidden)

    class Analysis:

        # The poseIdVec vector addresses an input poseid for each traj frame
        poseIdVec = OEField("PoseIdVec", Types.IntVec, meta=_metaHidden)

        # The OETraj Field is for the record containing Traj OEMols and energies
        oetraj_rec = OEField("OETraj", Types.Record, meta=_metaHidden)

        # The TrajIntE Field is for the record containing Traj interaction energies
        oeintE_rec = OEField("TrajIntE", Types.Record, meta=_metaHidden)

        # The TrajIntEDict Field is for the POD Dictionary containing Traj interaction energies
        oeintE_dict = OEField("TrajIntEDict", Types.JSONObject, meta=_metaHidden)

        # The TrajPBSA Field is for the record containing Traj PBSA energies
        oepbsa_rec = OEField("TrajPBSA", Types.Record, meta=_metaHidden)

        # The TrajPBSADict Field is for the POD Dictionary containing Traj PBSA energies
        oepbsa_dict = OEField("TrajPBSADict", Types.JSONObject, meta=_metaHidden)

        # The TrajClus Field is for the record containing Traj ligand clustering results
        oeclus_rec = OEField("TrajClus", Types.Record, meta=_metaHidden)

        # The TrajClusDict Field is for the POD Dictionary containing Traj ligand clustering results
        oeclus_dict = OEField("TrajClusDict", Types.JSONObject, meta=_metaHidden)

        # The ClusPopDict Field is for the POD Dictionary containing conf/cluster population results
        cluspop_dict = OEField("ClusPopDict", Types.JSONObject, meta=_metaHidden)

        # The AnalysesDone Field is for a list of the analyses that have been done
        analysesDone = OEField("AnalysesDone", Types.StringVec, meta=_metaHidden)

        # The Lig_Conf_Data Field is for the record containing Traj conf data for all confs
        oetrajconf_rec = OEField("Lig_Conf_Data", Types.RecordVec, meta=_metaHidden)

        # The vector of ligand Traj RMSDs from the initial pose
        lig_traj_rmsd = OEField('LigTrajRMSD', Types.FloatVec,
                                meta=OEFieldMeta().set_option(Meta.Units.Length.Ang))

        # The mmpbsa Field contains the vector of per-frame mmpbsa values over the whole trajectory
        zapMMPBSA_fld = OEField("OEZap_MMPBSA6_Bind", Types.FloatVec,
                                  meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))

        # mmpbsa ensemble average over the whole trajectory
        mmpbsa_traj_mean = OEField('MMPBSATrajMean', Types.Float,
                                   meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

        metaMMPBSA_traj_serr = OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol)
        metaMMPBSA_traj_serr.add_relation(Meta.Relations.ErrorsFor, mmpbsa_traj_mean)
        mmpbsa_traj_serr = OEField('MMPBSATrajSerr', Types.Float, meta=metaMMPBSA_traj_serr)

        # The number of major clusters found
        n_major_clusters = OEField("n major clusters", Types.Int)

        # Trajectory cluster averages and medians of protein and ligand
        ClusLigAvg_fld = OEField('ClusLigAvgMol', Types.Chem.Mol)
        ClusProtAvg_fld = OEField('ClusProtAvgMol', Types.Chem.Mol)
        ClusLigMed_fld = OEField('ClusLigMedMol', Types.Chem.Mol)
        ClusProtMed_fld = OEField('ClusProtMedMol', Types.Chem.Mol)

        max_waters = OEField("MaxWaters_OPLMD", Types.Int, meta=_metaHidden)

        # Free Energy
        # Analysis Fields
        free_energy = OEField('FE_OPLMD', Types.Float,
                              meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

        metaFreeEnergy_err = OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol)
        metaFreeEnergy_err.add_relation(Meta.Relations.ErrorsFor, free_energy)
        free_energy_err = OEField('FE_Error_OPLMD', Types.Float, meta=metaFreeEnergy_err)

    class FEC:
        # Free Energy
        free_energy = OEField('FE_OPLMD', Types.Float,
                              meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

        metaFreeEnergy_err = OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol)
        metaFreeEnergy_err.add_relation(Meta.Relations.ErrorsFor, free_energy)
        free_energy_err = OEField('FE_Error_OPLMD', Types.Float, meta=metaFreeEnergy_err)

        class RBFEC:
            # Oriented Edge field for relative free energy calculations
            # The first integer of the list is the ligand ID of the starting
            # thermodynamic state and the second the final one
            edgeid = OEField("EdgeID_OPLMD", Types.Int, meta=_metaHidden)
            edge_name = OEField("EdgeName_OPLMD", Types.String)

            # The Thermodynamics leg type is used for Bound and
            # UnBound State run identification
            thd_leg_type = OEField("Thd_Leg_OPLMD", Types.String, meta=_metaHidden)

            class NESC:

                state_A = OEField("StateA_OPLMD", Types.RecordVec, meta=_metaHidden)
                state_B = OEField("StateB_OPLMD", Types.RecordVec, meta=_metaHidden)

                forward = OEField("Forward_OPLMD", Types.String)
                reverse = OEField("Reverse_OPLMD", Types.String)

                direction = OEField("Direction_OPLMD", Types.String)

                gmx_top = OEField("GMX_Top_OPLMD", Types.String, meta=_metaHidden)
                gmx_gro = OEField("GMX_Gro_OPLMD", Types.String, meta=_metaHidden)

                work = OEField("GMX_Work_OPLMD", Types.Float,
                               meta=OEFieldMeta().set_option(Meta.Units.Energy.kJ_per_mol))
                frame_count = OEField("frame_count", Types.Int, meta=_metaHidden)

                # The Work record is used to collect the data related to the
                # Work Forward and Reverse for the Bound and Unbound States
                work_rec = OEField("Work_Record_OPLMD", Types.Record)

                # The Relative Binding Affinity record collects data for the
                # different analysis methods used to compute it
                DDG_rec = OEField("DDG_Record_OPLMD", Types.Record)

