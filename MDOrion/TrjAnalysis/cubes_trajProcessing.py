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

from orionplatform.mixins import RecordPortsMixin

from floe.api import (ParallelMixin,
                      parameters,
                      ComputeCube)

from snowball.utils.log_params import LogFieldParam

from MDOrion.Standards import Fields, MDStageNames

from oeommtools import utils as oeutils

from orionclient.session import in_orion

import MDOrion.TrjAnalysis.utils as utl

import oetrajanalysis.TrajMMPBSA_utils as mmpbsa

from MDOrion.TrjAnalysis.water_utils import nmax_waters

from openeye import oechem

import os,math

import traceback

from datarecord import (Types,
                        OEPrimaryMolField,
                        OEField,
                        OERecord)

from MDOrion.Standards.mdrecord import MDDataRecord

from MDOrion.Standards.standards import CollectionsNames

# use a really large float as a magic number to replace NaNs to avoid Orion WriterCube errors
magic_big_float_to_replace_NaN = 4.0e+256


class TrajToOEMolCube(RecordPortsMixin, ComputeCube):
    title = 'Traj to OEMol Cube'
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Trajectory', 'Ligand', 'Protein']

    description = """
    Converting MD Traj into multiconf OEMols for Ligand and Protein.
    This Cube will take in the MD traj file containing
    the solvated protein:ligand complex and extract
    multiconf OEMols for Ligand and Protein.
    """

    uuid = "3ad0e991-712f-4a87-903e-4e0edc774bb3"

    # for Exception Handler
    log_field = LogFieldParam()

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 32000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    water_cutoff = parameters.DecimalParameter(
        'water_cutoff',
        default=15.0,
        help_text="""The cutoff distance in angstroms to select waters around the
        protein-ligand binding site for each trajectory frame""")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        return

    def process_failed(self, record, port, last_error):
        print("Failed to process record", record, flush=True)
        self.failure.emit(record)

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Initialize system_title for exception handling
            system_title = ''

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            # Logger string
            opt['Logger'].info(' ')
            system_title = mdrecord.get_title
            opt['Logger'].info('{}: Attempting MD Traj conversion into OEMols'.format(system_title))

            traj_fn = mdrecord.get_stage_trajectory()

            # opt['Logger'].info('{} Temp Directory: {}'.format(system_title, os.path.dirname(traj_fn)))
            # opt['Logger'].info('{} Trajectory filename: {}'.format(system_title, traj_fn))

            # Generate multi-conformer protein and ligand OEMols from the trajectory
            opt['Logger'].info('{} Generating trajectory OEMols'.format(system_title))

            flask = mdrecord.get_flask

            md_components = mdrecord.get_md_components

            # Checking
            if not md_components.has_ligand:
                raise ValueError("The Ligand MD component is missing and the Analysis cannot be performed")
            if not md_components.has_protein:
                raise ValueError("The Protein MD component is missing and the Analysis cannot be performed")
            if not md_components.has_water:
                raise ValueError("The Water MD component is missing and the Analysis cannot be performed")

            if not record.has_value(Fields.ligand):
                raise ValueError("The Reference Ligand cannot be found. This Cube should be used "
                                 "after the Solvate and Run Protein and Ligand MD floe has been completed")

            # Check Ligand Isomeric Smiles
            lig_ref = record.get_value(Fields.ligand)
            lig_comp = md_components.get_ligand

            smi_lig_comp = oechem.OECreateSmiString(lig_comp)
            smi_lig_ref = oechem.OECreateSmiString(lig_ref)

            if smi_lig_ref != smi_lig_comp:
                raise ValueError("The Reference Ligand and the Ligand MD component mismatch."
                                 " Isomeric Smiles String check failure: {} vs {}".format(smi_lig_comp, smi_lig_ref))

            ptraj, ltraj, wtraj = utl.extract_aligned_prot_lig_wat_traj(md_components, flask, traj_fn, opt,
                                                                        water_cutoff=opt['water_cutoff'])

            if not record.has_value(Fields.ligand_name):
                raise ValueError("The ligand name is missing from the record")

            if not record.has_value(Fields.protein_name):
                raise ValueError("The protein name is missing from the record")

            ltraj.SetTitle(record.get_value(Fields.ligand_name))
            ptraj.SetTitle(record.get_value(Fields.protein_name))

            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'.format(
                system_title, ptraj.NumAtoms(), ptraj.NumConfs()))
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'.format(
                system_title, ltraj.NumAtoms(), ltraj.NumConfs()))
            opt['Logger'].info('{} #atoms, #confs in water traj OEMol: {}, {}'.format(
                system_title, wtraj.NumAtoms(), wtraj.NumConfs()))

            # Create new record with OETraj results
            oetrajRecord = OERecord()

            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ltraj)

            if wtraj:
                oetrajRecord.set_value(OEField('WatTraj', Types.Chem.Mol), wtraj)

            if in_orion():
                oetrajRecord.set_value(Fields.collections, {CollectionsNames.md: mdrecord.collection_id})

            mdrecord_traj = MDDataRecord(oetrajRecord)

            mdrecord_traj.set_protein_traj(ptraj, shard_name="ProteinTrajConfs_")

            record.set_value(Fields.Analysis.oetraj_rec, oetrajRecord)

            # update or initiate the list of analyses that have been done
            if record.has_value(Fields.Analysis.analysesDone):
                analysesDone = utl.RequestOEFieldType(record,  Fields.Analysis.analysesDone)
                analysesDone.append('OETraj')
            else:
                analysesDone = ['OETraj']

            record.set_value(Fields.Analysis.analysesDone, analysesDone)

            opt['Logger'].info('{}: saved protein, ligand  and water traj OEMols'.format(system_title))

            self.success.emit(record)

            del mdrecord
            del mdrecord_traj

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.log.error(traceback.format_exc())
            # Write field for Exception Handler
            msg = '{}: {} Cube: {}'.format(system_title, self.title, str(e))
            record.set_value(self.args.log_field, msg)
            # Return failed mol
            self.failure.emit(record)

        return


class ConcatenateTrajMMPBSACube(RecordPortsMixin, ComputeCube):
    title = "Concatenate Trajectory MMPBSA Energy Components"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['OEChem', 'TrajAnalysis', 'MMPBSA']
    description = """
    Protein-ligand trajectory MMPBSA interaction energies are concatenated based on
    pre-existing per-conformer vectors.

    Input:
    -------
    Data Record with the per-conformer trajectory MMPBSA interaction energy components.

    Output:
    -------
    Data Record - The various energy components associated with the Poisson-Boltzmann and
    Surface Area energies are attached to the record as per-frame vectors of floats.
    The energy units are in kcal/mol.
    """

    uuid = "4a36f62f-9eb2-4589-b884-aa4ca653c5cd"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            # Logger string
            opt['Logger'].info(' Beginning ConcatenateTrajMMPBSACube')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to concatenate per-conf MMPBSA energies'.format(system_title))

            # Go find the conformer records
            if not record.has_field(Fields.Analysis.oetrajconf_rec):
                raise ValueError('{} could not find the conformer record'.format(system_title))
            else:
                opt['Logger'].info('{} found the conformer record'.format(system_title))
            list_conf_rec = record.get_value(Fields.Analysis.oetrajconf_rec)

            # combine all data from the conf traj PBSA data for all conformers
            PBSAdata = dict()
            for confrec in list_conf_rec:
                confid = utl.RequestOEFieldType(confrec, Fields.confid)
                if not confrec.has_field(Fields.Analysis.oepbsa_dict):
                    raise ValueError('{} could not find the conf traj PBSA data for confid {}'.
                                     format(system_title, confid))
                confPBSAdata = confrec.get_value(Fields.Analysis.oepbsa_dict)

                for key in confPBSAdata.keys():
                    if key not in PBSAdata.keys():
                        PBSAdata[key] = confPBSAdata[key]
                    else:
                        PBSAdata[key] += confPBSAdata[key]

            #for key in PBSAdata.keys():
            #    opt['Logger'].info('ConcatenateTrajMMPBSACube PBSAdata[{}] length {}'.format(key, len(PBSAdata[key])))

            # Add the PBSAdata dict to the parent record
            record.set_value(Fields.Analysis.oepbsa_dict, PBSAdata)

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} in TrajPBSACube'.format(str(e)))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class TrajPBSACube(RecordPortsMixin, ComputeCube):
    title = "Trajectory Poisson-Boltzmann and Surface Area Energies"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['OEChem', 'Zap', 'TrajAnalysis', 'MMPBSA']
    description = """
    Protein-ligand interaction solvation energies are calculated on an existing MD trajectory.
    The trajectory is taken from pre-existing protein and ligand trajectory OEMols.
    The Poisson-Boltzmann and Surface Area methods in the OEZap toolkits are  used.
    The various energy components associated with the Poisson-Boltzmann and
    Surface Area energies are attached to the record as per-frame vectors of floats.
    The energy units are in kcal/mol.
    """

    uuid = "f6c96295-51fd-42df-8763-0f3b6f6d0e0d"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    explicit_water = parameters.BooleanParameter(
        'explicit_water',
        default=False,
        help_text="""Enable MMPBSA calculation with explicit water""")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            # Logger string
            opt['Logger'].info(' Beginning TrajPBSACube')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to compute MD Traj PBSA energies'.format(system_title))

            # Check that the OETraj analysis has been done
            analysesDone = utl.RequestOEFieldType(record, Fields.Analysis.analysesDone)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title))
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title))

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = utl.RequestOEFieldType(record, Fields.Analysis.oetraj_rec)
            opt['Logger'].info('{} found OETraj record'.format(system_title))
            ligTraj = utl.RequestOEField(oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                               .format(system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()))

            mdtrajrecord = MDDataRecord(oetrajRecord)

            if self.opt['explicit_water']:

                water_traj = oetrajRecord.get_value(OEField('WatTraj', Types.Chem.Mol))
                opt['Logger'].info('{} #atoms, #confs in water traj OEMol: {}, {}'
                                   .format(system_title, water_traj.NumAtoms(), water_traj.NumConfs()))

                protTraj = mdtrajrecord.get_protein_traj

                prot_wat = oechem.OEMol(protTraj.GetActive())
                oechem.OEAddMols(prot_wat, water_traj.GetActive())

                prot_wat.DeleteConfs()

                for pr_conf, wat_conf in zip(protTraj.GetConfs(), water_traj.GetConfs()):
                    pr_wat_conf = oechem.OEMol(pr_conf)
                    oechem.OEAddMols(pr_wat_conf, wat_conf)
                    pr_wat_conf_xyz = oechem.OEFloatArray(prot_wat.NumAtoms() * 3)
                    pr_wat_conf.GetCoords(pr_wat_conf_xyz)
                    prot_wat.NewConf(pr_wat_conf_xyz)

                protTraj = prot_wat
            else:
                protTraj = mdtrajrecord.get_protein_traj

            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                               .format(system_title, protTraj.NumAtoms(), protTraj.NumConfs()))

            # Compute PBSA energies for the protein-ligand complex
            PBSAdata = mmpbsa.TrajPBSA(ligTraj, protTraj)
            if PBSAdata is None:
                raise ValueError('{} Calculation of PBSA energies failed'.format(system_title))

            # generate Surface Areas energy for buried SA based on 0.006 kcal/mol/A^2
            PBSAdata['OEZap_SA6_Bind'] = [sa * -0.006 for sa in PBSAdata['OEZap_BuriedArea']]

            # If the OETraj Interaction Energies has been done calculate MMPBSA values
            if 'TrajIntE' in analysesDone:
                opt['Logger'].info('{} found TrajIntE analyses'.format(system_title) )

                # Extract the relevant P-L Interaction Energies from the record
                intEdata = record.get_value(Fields.Analysis.oeintE_dict)
                opt['Logger'].info('{} found Traj intEdata data'.format(system_title))

                if self.opt['explicit_water']:

                    PLIntE = intEdata['protein_and_water_ligand_interE']
                    opt['Logger'].info('{} found Protein-Water and Ligand force field interaction energies'
                                       .format(system_title))
                else:

                    PLIntE = intEdata['protein_ligand_interE']
                    opt['Logger'].info('{} found Protein-Ligand force field interaction energies'
                                       .format(system_title))

                # Calculate  and store MMPB and MMPBSA energies on the trajPBSA record
                PBSAdata['OEZap_MMPB_Bind'] = [eInt+eDesol for eInt, eDesol in
                                               zip(PLIntE, PBSAdata['OEZap_PB_Desolvation'])]
                PBSAdata['OEZap_MMPBSA6_Bind'] = [eMMPB+eSA6 for eMMPB,eSA6 in
                                                  zip(PBSAdata['OEZap_MMPB_Bind'], PBSAdata['OEZap_SA6_Bind'])]

            # list field and change any NaNs to a really big float
            for key in PBSAdata.keys():
                opt['Logger'].info('{} TrajPBSACube PBSAdata[{}] of length {}'
                                   .format(system_title,key,len(PBSAdata[key])) )
                # change any NaNs to a really big float or else Orion WriterCube fails on JSON dict
                for i, x in enumerate(PBSAdata[key]):
                    if math.isnan(x):
                        opt['Logger'].info('{} found a NaN at PBSAdata[{}][{}]'.format(system_title,key,i))
                        PBSAdata[key][i] = magic_big_float_to_replace_NaN

            # Add the PBSAdata dict to the record
            record.set_value(Fields.Analysis.oepbsa_dict, PBSAdata)

            analysesDone.append('TrajPBSA')
            record.set_value(Fields.Analysis.analysesDone, analysesDone)
            opt['Logger'].info('{} finished writing TrajPBSA OERecord'.format(system_title))

            self.success.emit(record)

            del mdtrajrecord

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} in TrajPBSACube'.format(str(e)))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class TrajInteractionEnergyCube(RecordPortsMixin, ComputeCube):
    title = "Trajectory Interaction Energies"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['OEChem', 'OpenMM', 'TrajAnalysis', 'MMPBSA']
    description = """
    Protein-ligand interaction energies are calculated on an existing MD trajectory.
    The trajectory is taken from pre-existing protein and ligand trajectory OEMols.
    The forcefield used is taken from the parmed object associated with the trajectory
    OEMols. The various energy components associated with the protein-ligand
    interaction energies are attached to the record as per-frame vectors of floats.
    This includes the MM interaction potential energies and their ligand, protein
    and complex components. The energy units are in kcal/mol.
    """

    uuid = "d10a770d-fcd2-4d09-bf00-00d6a00353de"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 32000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            # Logger string
            opt['Logger'].info(' Beginning TrajInteractionEnergyCube')

            mdrecord = MDDataRecord(record)

            system_title = mdrecord.get_title

            opt['Logger'].info('{} Attempting to compute MD Traj protein-ligand Interaction energies'
                               .format(system_title))

            # Check that the OETraj analysis has been done
            analysesDone = utl.RequestOEFieldType(record, Fields.Analysis.analysesDone)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title))
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title))

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = utl.RequestOEFieldType(record, Fields.Analysis.oetraj_rec)
            opt['Logger'].info('{} found OETraj record'.format(system_title))
            ligTraj = utl.RequestOEField(oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                               .format(system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()))

            mdtrajrecord = MDDataRecord(oetrajRecord)
            protTraj = mdtrajrecord.get_protein_traj

            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'.
                               format(system_title, protTraj.NumAtoms(), protTraj.NumConfs()))

            water_traj = oetrajRecord.get_value(OEField('WatTraj', Types.Chem.Mol))
            opt['Logger'].info('{} #atoms, #confs in water traj OEMol: {}, {}'
                               .format(system_title, water_traj.NumAtoms(), water_traj.NumConfs()))

            prmed = mdrecord.get_parmed(sync_stage_name='last')

            # Compute interaction energies for the protein, ligand, complex and water subsystems
            intEdata = mmpbsa.ProtLigWatInteractionEFromParmedOETraj(prmed, ligTraj, protTraj, water_traj, opt)

            if intEdata is None:
                raise ValueError('{} Calculation of Interaction Energies failed'.format(system_title))

            # protein and ligand traj OEMols now have parmed charges on them; save these
            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ligTraj)
            record.set_value(Fields.Analysis.oetraj_rec, oetrajRecord)

            # list the energy terms in the intEdata dict to be stored on the record
            for key in intEdata.keys():
                opt['Logger'].info('{} traj intEdata[{}] of length {}'
                                   .format(system_title,key,len(intEdata[key])) )
                # change any NaNs to a really big float or else Orion WriterCube fails on JSON dict
                for i, x in enumerate(intEdata[key]):
                    if math.isnan(x):
                        opt['Logger'].info('{} found a NaN at intEdata[{}][{}]'.format(system_title,key,i))
                        intEdata[key][i] = magic_big_float_to_replace_NaN

            # Add the intEdata dict to the record
            record.set_value(Fields.Analysis.oeintE_dict, intEdata)

            # Add the trajIntE record to the parent record
            #record.set_value(Fields.Analysis.oeintE_rec, trajIntE)

            analysesDone.append('TrajIntE')
            record.set_value(Fields.Analysis.analysesDone, analysesDone)
            opt['Logger'].info('{} finished writing trajIntE OERecord'.format(system_title) )

            self.success.emit(record)

            del mdrecord
            del mdtrajrecord

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in TrajInteractionEnergyCube on {}'.format(str(e),system_title) )
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class ConfTrajsToLigTraj(RecordPortsMixin, ComputeCube):
    title = 'Conf Trajs To Ligand Traj'
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Clustering', 'Ligand', 'Protein']

    description = """
    Combine individual conformer trajectory OEMols into single ligand OEMol

    This Cube will read in and combine the individual MD traj OEMols for each conformer
    into a single MD traj OEMol for the whole ligand, in preparation for clistering.
    It will do this for both the protein and ligand components of the complex.
    """

    uuid = "16e88f36-28fe-4de2-8e61-bb0e5efff96f"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Logger string
            opt['Logger'].info(' Beginning ConfTrajsToLigTraj')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to combine conf traj OEMols into ligand traj OEMol'
                               .format(system_title))

            # Go find the ligand and LigTraj fields in each of the conformer records
            if not record.has_field(Fields.Analysis.oetrajconf_rec):
                raise ValueError('{} could not find the conformer record'.format(system_title))
            else:
                opt['Logger'].info('{} found the conformer record'.format(system_title))

            # set up ligand and LigTraj lists then loop over conformer records
            poseIdVec = []
            ligTrajConfs = []
            protTrajConfs = []
            watTrajConfs = []
            list_conf_rec = record.get_value(Fields.Analysis.oetrajconf_rec)
            for confrec in list_conf_rec:
                confid = utl.RequestOEFieldType(confrec, Fields.confid)

                if not confrec.has_field(Fields.Analysis.oetraj_rec):
                    raise ValueError('{} confID {}: could not find traj record'.format(system_title,confid))
                oetrajRecord = confrec.get_value(Fields.Analysis.oetraj_rec)

                # Extract the ligand traj OEMol from the OETraj record
                ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
                poseIdVec += [confid]*ligTraj.NumConfs()
                ligTrajConfs.append(ligTraj)
                opt['Logger'].info('{} confID {}: adding ligTraj with {} atoms, {} confs'.format(
                    system_title, confid, ligTraj.NumAtoms(), ligTraj.NumConfs()) )

                # Extract the activeSite water traj OEMol from the OETraj record
                watTraj = utl.RequestOEField( oetrajRecord, 'WatTraj', Types.Chem.Mol)
                watTrajConfs.append(watTraj)
                opt['Logger'].info('{} confID {}: adding watTraj with {} atoms, {} confs'.format(
                    system_title, confid, watTraj.NumAtoms(), watTraj.NumConfs()) )

                # Extract the protTraj OEMol from the OETraj record
                mdtrajrecord = MDDataRecord(oetrajRecord)
                protTraj = mdtrajrecord.get_protein_traj
                protTrajConfs.append(protTraj)
                opt['Logger'].info('{} confID {}: adding protTraj with {} atoms, {} confs'.format(
                    system_title, confid, protTraj.NumAtoms(), protTraj.NumConfs()) )
                del mdtrajrecord

            if len(ligTrajConfs) < 1 or len(protTrajConfs) < 1:
                raise ValueError('{} empty list of lig or protein trajectory OEMols'.format(system_title))

            ligTraj = oechem.OEMol(ligTrajConfs[0])
            xyz = oechem.OEFloatArray(3*ligTraj.GetMaxAtomIdx())
            for trajMol in ligTrajConfs[1:]:
                for conf in trajMol.GetConfs():
                    conf.GetCoords(xyz)
                    ligTraj.NewConf(xyz)
            opt['Logger'].info('{} composite ligTraj has {} atoms, {} confs'.format(
                system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )

            watTraj = oechem.OEMol(watTrajConfs[0])
            xyz = oechem.OEFloatArray(3*watTraj.GetMaxAtomIdx())
            for trajMol in watTrajConfs[1:]:
                for conf in trajMol.GetConfs():
                    conf.GetCoords(xyz)
                    watTraj.NewConf(xyz)
            opt['Logger'].info('{} composite watTraj has {} atoms, {} confs'.format(
                system_title, watTraj.NumAtoms(), watTraj.NumConfs()) )

            protTraj = protTrajConfs[0]
            xyz = oechem.OEFloatArray(3*protTraj.GetMaxAtomIdx())
            for trajMol in protTrajConfs[1:]:
                for conf in trajMol.GetConfs():
                    conf.GetCoords(xyz)
                    protTraj.NewConf(xyz)
            opt['Logger'].info('{} composite protTraj has {} atoms, {} confs'.format(
                system_title, protTraj.NumAtoms(), protTraj.NumConfs()))

            record.set_value(Fields.Analysis.poseIdVec, poseIdVec)

            # Create new record with OETraj results
            oetrajRecord = OERecord()
            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ligTraj)
            if watTraj:
                oetrajRecord.set_value(OEField('WatTraj', Types.Chem.Mol), watTraj)

            if in_orion():
                collection_dic = utl.RequestOEFieldType(record, Fields.collections)
                oetrajRecord.set_value(Fields.collections, collection_dic)
            mdrecord_traj = MDDataRecord(oetrajRecord)
            mdrecord_traj.set_protein_traj(protTraj, shard_name="ProteinTrajConfs_")

            record.set_value(Fields.Analysis.oetraj_rec, oetrajRecord)

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ConfTrajsToLigTraj on {}'.format(str(e), system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class ConformerGatheringData(RecordPortsMixin, ComputeCube):
    title = "MD Conformer Gathering Data"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Ligand', 'Protein']

    description = """
    This Cube gathers together conformers related to the same ligand and their information
    in a new record containing the multi conformer ligand and each conformer record info.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 32000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    uuid = "1cc827a6-a2c2-4b51-a705-dd082f3a200c"

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

        # This dictionary stores for each ligand all its conformers in a list
        self.lig_sys_ids = dict()

    def process(self, record, port):
        try:
            mdrecord = MDDataRecord(record)

            sys_id = mdrecord.get_lig_id

            if sys_id not in self.lig_sys_ids.keys():
                self.lig_sys_ids[sys_id] = [record]
            else:
                self.lig_sys_ids[sys_id].append(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):

        try:
            for sys_id, list_conf_rec in self.lig_sys_ids.items():

                # Save the first record to emit in failure cases
                self.record = list_conf_rec

                # catch case where for some reason the conf list list_conf_rec is empty
                if len(list_conf_rec) < 1:
                    print('{} does not have any conformer data'.format(sys_id) )
                    continue
                elif len(list_conf_rec) > 1:
                    # Conformers for each ligand are sorted based on their confid in each ligand record
                    list_conf_rec.sort(key=lambda x: x.get_value(Fields.confid))

                new_rec = OERecord()
                new_rec.set_value(Fields.Analysis.oetrajconf_rec, list_conf_rec)
                # Get the first conf to move some general ligand data up to the top level
                rec0 = list_conf_rec[0]
                #   copy all the initial fields in Fields.ligInit_rec up to the top level
                init_rec = rec0.get_value(Fields.ligInit_rec)

                # TODO METADATA IS NOT COPIED?
                for field in init_rec.get_fields():
                    new_rec.set_value(field, init_rec.get_value(field))
                #   next, fields that will simply be copied and not further used here
                protein = rec0.get_value(Fields.protein)
                new_rec.set_value(Fields.protein, protein)
                ligid = rec0.get_value(Fields.ligid)
                new_rec.set_value(Fields.ligid, ligid)
                if in_orion():
                    collection_dic = rec0.get_value(Fields.collections)
                    new_rec.set_value(Fields.collections, collection_dic)
                #   finally, fields that will be copied and also further used here
                lig_multi_conf = oechem.OEMol(rec0.get_value(Fields.ligand))
                protein_name = rec0.get_value(Fields.protein_name)

                # MD Components copied at the ligand top level
                new_rec.set_value(Fields.md_components, rec0.get_value(Fields.md_components))

                # if >1 confs, add their confs to the parent ligand at the top level
                for rec in list_conf_rec[1:]:
                    lig_multi_conf.NewConf(rec.get_value(Fields.ligand))

                # get name of initial molecule
                if new_rec.has_value(OEPrimaryMolField()):
                    init_mol = new_rec.get_value(OEPrimaryMolField())
                else:
                    print('{} ConformerGatheringData: new_rec cannot find the OEPrimaryMolField'.format(sys_id) )
                    continue
                lig_title = init_mol.GetTitle()
                lig_multi_conf.SetTitle(lig_title)
                # regenerate protein-ligand title since all titles on conformers include conformer id
                title = 'p' + protein_name + '_l' + lig_title
                # set other fields on the new record
                new_rec.set_value(Fields.title, title)
                new_rec.set_value(Fields.ligand, lig_multi_conf)
                new_rec.set_value(Fields.primary_molecule, lig_multi_conf)
                new_rec.set_value(Fields.protein_name, protein_name)
                new_rec.set_value(Fields.ligand_name, lig_title)

                self.success.emit(new_rec)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(self.record)


class NMaxWatersLigProt(RecordPortsMixin, ComputeCube):
    title = "NMax Waters"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Ligand', 'Protein', 'Waters']

    description = """
    This Cube determines the max number of waters for all the ligands that
    fits between the protein and ligand molecular surfaces. The cutoff distance
    parameters determines the max distance used between the volume grid points
    and the ligand-protein.
    """

    uuid = "c8608012-dc09-47de-8d23-fe2338419ff3"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    cutoff = parameters.DecimalParameter(
        'cutoff',
        default=5.0,
        help_text="Cutoff Distance between Volume grid points and ligand-protein in A")

    explicit_water = parameters.BooleanParameter(
        'explicit_water',
        default=False,
        help_text="""Enable MMPBSA calculation with explicit water""")

    water_number = parameters.IntegerParameter(
        'water_number',
        default=0,
        help_text="""If different from zero the selected water number will be used"""
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.nwaters = list()
        self.records = list()

    def process(self, record, port):
        try:
            mdrecord = MDDataRecord(record)

            protein = mdrecord.get_protein

            protein, ligand, water, exc = oeutils.split(protein, ligand_res_name='LIG')

            if protein.NumAtoms() == 0:
                raise ValueError("The Protein Atom number is zero")

            ligand = mdrecord.get_ligand

            if self.opt['explicit_water']:

                if self.opt['water_number'] != 0:
                    self.opt['Logger'].info("User selected number of waters: {}".format(self.opt['water_number']))
                    nmax = self.opt['water_number']
                else:
                    nmax = nmax_waters(protein, ligand, self.opt['cutoff'])
            else:
                self.opt['Logger'].info("MMPBSA Explicit Water set off")
                nmax = 0

            self.nwaters.append(nmax)
            self.records.append(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):

        max_waters = max(self.nwaters)

        if max_waters == 0:
            self.opt['Logger'].warn("[{}] Max number of waters is zero".format(self.title))
        else:
            self.opt['Logger'].info("[{}] Max number of Waters: {}".format(self.title, max_waters))

        for rec in self.records:
            rec.set_value(Fields.Analysis.max_waters, max_waters)
            self.success.emit(rec)

        return


class ParallelTrajToOEMolCube(ParallelMixin, TrajToOEMolCube):
    title = "Parallel " + TrajToOEMolCube.title
    description = "(Parallel) " + TrajToOEMolCube.description
    uuid = "c6435bfa-9e2c-4c6d-a59f-81baff1dc7b8"


class ParallelTrajInteractionEnergyCube(ParallelMixin, TrajInteractionEnergyCube):
    title = "Parallel " + TrajInteractionEnergyCube.title
    description = "(Parallel) " + TrajInteractionEnergyCube.description
    uuid = "a6a11dbb-bc25-4548-bf1a-471bda2f0406"


class ParallelConfTrajsToLigTraj(ParallelMixin, ConfTrajsToLigTraj):
    title = "Parallel " + ConfTrajsToLigTraj.title
    description = "(Parallel) " + ConfTrajsToLigTraj.description
    uuid = "5b81622a-fbc9-47e1-a8f2-6a7e5675df8a"


class ParallelConcatenateTrajMMPBSACube(ParallelMixin, ConcatenateTrajMMPBSACube):
    title = "Parallel " + ConcatenateTrajMMPBSACube.title
    description = "(Parallel) " + ConcatenateTrajMMPBSACube.description
    uuid = "ca6b8b5d-a7d1-494b-ae12-253ea96ecabf"


class ParallelTrajPBSACube(ParallelMixin, TrajPBSACube):
    title = "Parallel " + TrajPBSACube.title
    description = "(Parallel) " + TrajPBSACube.description
    uuid = "1f62c9ee-1d83-469c-b6a0-3b20dc64a1e1"


