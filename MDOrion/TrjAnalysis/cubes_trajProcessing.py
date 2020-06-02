from orionplatform.mixins import RecordPortsMixin

from floe.api import (ParallelMixin,
                      parameters,
                      ComputeCube)

from MDOrion.Standards import Fields, MDStageNames

from oeommtools import utils as oeutils

from floereport import FloeReport, LocalFloeReport

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

from os import environ

import MDOrion.TrjAnalysis.utils as utl

import MDOrion.TrjAnalysis.TrajMMPBSA_utils as mmpbsa

from MDOrion.TrjAnalysis.water_utils import nmax_waters

import oetrajanalysis.OETrajBasicAnalysis_utils as oetrjutl

import ensemble2img

from tempfile import TemporaryDirectory

from openeye import oechem

import oetrajanalysis.Clustering_utils as clusutl

from openeye import oedepict

import os

import traceback

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from MDOrion.Standards.mdrecord import MDDataRecord

from MDOrion.TrjAnalysis.TrajAnFloeReport_utils import (_clus_floe_report_header,
                                                        _clus_floe_report_header2,
                                                        _clus_floe_report_midHtml0,
                                                        _clus_floe_report_midHtml1,
                                                        _clus_floe_report_midHtml2,
                                                        _clus_floe_report_stripPlots,
                                                        _clus_floe_report_Trailer,
                                                        trim_svg,
                                                        MakeClusterInfoText)


class TrajToOEMolCube(RecordPortsMixin, ComputeCube):
    title = 'Traj to OEMol Cube'
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Trajectory', 'Ligand', 'Protein']

    description = """
    Converting MD Traj into multiconf OEMols for Ligand and Protein.
    This cube will take in the MD traj file containing
    the solvated protein:ligand complex and extract
    multiconf OEMols for Ligand and Protein.
    """

    uuid = "3ad0e991-712f-4a87-903e-4e0edc774bb3"

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

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            # Logger string
            opt['Logger'].info(' ')
            system_title = mdrecord.get_title
            sys_id = mdrecord.get_flask_id
            opt['Logger'].info('{}: Attempting MD Traj conversion into OEMols'.format(system_title))

            traj_fn = mdrecord.get_stage_trajectory()

            opt['Logger'].info('{} Temp Directory: {}'.format(system_title, os.path.dirname(traj_fn)))
            opt['Logger'].info('{} Trajectory filename: {}'.format(system_title, traj_fn))

            setupOEMol = mdrecord.get_stage_topology(stg_name=MDStageNames.ForceField)

            opt['Logger'].info('{} Setup topology has {} atoms'.format(system_title, setupOEMol.NumAtoms()))

            # Generate multi-conformer protein and ligand OEMols from the trajectory
            opt['Logger'].info('{} Generating protein and ligand trajectory OEMols'.format(system_title))

            flask = mdrecord.get_flask

            ptraj, ltraj, wtraj = utl.extract_aligned_prot_lig_wat_traj(setupOEMol, flask, traj_fn, opt,
                                                                        water_cutoff=opt['water_cutoff'])
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
                oetrajRecord.set_value(Fields.collection, mdrecord.collection_id)

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

    # uuid = "f6c96295-51fd-42df-8763-0f3b6f6d0e0d"

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

            # find the PBSAdata dict from inside the first conf
            #if not rec0.has_field(Fields.Analysis.oepbsa_dict):
            #    raise ValueError('{} could not find the PBSAdata dict for confID {}'.format(system_title, confid))
            #opt['Logger'].info('{}: found PBSAdata dict for confID {}'.format(system_title, confid))
            #PBSAdata = rec0.get_value(Fields.Analysis.oepbsa_dict)
            #for key in PBSAdata:
            #    opt['Logger'].info('{}: PBSAdata {} : {} values'.format(system_title, key, len(PBSAdata[key])))

            # make a dict of all FloatVec fields from inside the trajPBSA_rec record of the first conf
            PBSAdata = dict()
            for confrec in list_conf_rec:
                confid = utl.RequestOEFieldType(confrec, Fields.confid)
                if not confrec.has_field(Fields.Analysis.oepbsa_rec):
                    raise ValueError('{} could not find the trajPBSA record for confid {}'.
                                     format(system_title, confid))
                trajPBSA_rec = confrec.get_value(Fields.Analysis.oepbsa_rec)
                for field in trajPBSA_rec.get_fields():
                    fname = field.get_name()
                    if field.get_type() is Types.FloatVec:
                        opt['Logger'].info(
                            'ConcatenateTrajMMPBSACube adding FloatVec {}'.format(fname))
                        if fname not in PBSAdata.keys():
                            PBSAdata[fname] = trajPBSA_rec.get_value(field)
                        else:
                            PBSAdata[fname] += trajPBSA_rec.get_value(field)
                    else:
                        opt['Logger'].info('ConcatenateTrajMMPBSACube skipping field {}'.format(fname))
            for key in PBSAdata.keys():
                opt['Logger'].info('ConcatenateTrajMMPBSACube PBSAdata[{}] length {}'.format(key, len(PBSAdata[key])))

            # Add the trajPBSA record to the parent record
            record.set_value(Fields.Analysis.oepbsa_dict, PBSAdata)



            # find the trajPBSA record from inside the first conf
            rec0 = list_conf_rec[0]
            confid = utl.RequestOEFieldType(rec0, Fields.confid)
            if not rec0.has_field(Fields.Analysis.oepbsa_rec):
                raise ValueError('{} could not find the trajPBSA record for confID {}'.format(system_title, confid))
            trajPBSA_rec = rec0.get_value(Fields.Analysis.oepbsa_rec)
            opt['Logger'].info('{}: found trajPBSA record for confID {}'.format(system_title, confid))

            # make a list of all FloatVec fields from inside the trajPBSA_rec record of the first conf
            floatVecFields = []
            for field in trajPBSA_rec.get_fields():
                if field.get_type() is Types.FloatVec:
                    opt['Logger'].info('ConcatenateTrajMMPBSACube adding FloatVec field {}'.format(field.get_name()))
                    floatVecFields.append(field)
                else:
                    opt['Logger'].info('ConcatenateTrajMMPBSACube skipping field {}'.
                                       format(field.get_name()))
            if len(floatVecFields)<1:
                raise ValueError('{} ConcatenateTrajMMPBSACube found no FloatVec Fields'.format(system_title))

            # for each field, concatenate the floatvecs from all confs and put in new record
            new_rec = OERecord()
            for field in floatVecFields:
                vec = []
                for confrec in list_conf_rec:
                    confid = utl.RequestOEFieldType(confrec, Fields.confid)
                    if not confrec.has_field(Fields.Analysis.oepbsa_rec):
                        raise ValueError('{} could not find the trajPBSA record for confid {}'.
                                         format(system_title,confid))
                    trajPBSA_rec = confrec.get_value(Fields.Analysis.oepbsa_rec)
                    vec += trajPBSA_rec.get_value(field)
                new_rec.set_value(field, vec)
                #opt['Logger'].info('ConcatenateTrajMMPBSACube added FloatVec {} of size {}'.format(field.get_name(),len(vec)))

            # Add the trajPBSA record to the parent record
            record.set_value(Fields.Analysis.oepbsa_rec, new_rec)

            # Get concatenated MMPBSA values for overall mean and std err
            zapMMPBSA = new_rec.get_value(Fields.Analysis.zapMMPBSA_fld, vec)
            # Clean average MMPBSA to avoid nans and high zap energy values
            avg_mmpbsa, serr_mmpbsa = utl.clean_mean_serr(zapMMPBSA)

            # Add to the record the MMPBSA mean and std
            record.set_value(Fields.Analysis.mmpbsa_traj_mean, avg_mmpbsa)
            record.set_value(Fields.Analysis.mmpbsa_traj_serr, serr_mmpbsa)
            # Add to the record the Average MMPBSA floe report label
            record.set_value(Fields.floe_report_label, "MMPBSA score:<br>{:.1f}  &plusmn; {:.1f} kcal/mol".
                             format(avg_mmpbsa, serr_mmpbsa))

            opt['Logger'].info('{} finished ConcatenateTrajMMPBSACube'.format(system_title))

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
            zapBind, zapBindPB, zapDesolEl, zapIntEl, zapBindSA25, saBuried = mmpbsa.TrajPBSA(
                                       ligTraj, protTraj)
            if zapBind is None:
                raise ValueError('{} Calculation of PBSA energies failed'.format(system_title))

            # generate Surface Areas energy for buried SA based on 0.006 kcal/mol/A^2
            zapBindSA6 = [sa * -0.006 for sa in saBuried]

            # make a dict of these energy terms to store on the record
            PBSAdata = dict()
            PBSAdata['OEZap_PBSA25_Bind'] = zapBind
            PBSAdata['OEZap_PB_Bind'] = zapBindPB
            PBSAdata['OEZap_PB_Desolvation'] = zapDesolEl
            PBSAdata['OEZap_PB_Interaction'] = zapIntEl
            PBSAdata['OEZap_SA25_Bind'] = zapBindSA25
            PBSAdata['OEZap_BuriedArea'] = saBuried
            PBSAdata['OEZap_SA6_Bind'] = zapBindSA6

            # Create new record with traj interaction energy results
            opt['Logger'].info('{} writing trajPBSA OERecord'.format(system_title) )
            trajPBSA = OERecord()

            zapBind_field = OEField("OEZap_PBSA25_Bind", Types.FloatVec,
                                    meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBind_field, zapBind)

            zapBindPB_field = OEField("OEZap_PB_Bind", Types.FloatVec,
                                      meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindPB_field, zapBindPB)

            zapDesolEl_field = OEField("OEZap_PB_Desolvation", Types.FloatVec,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapDesolEl_field, zapDesolEl)

            zapIntEl_field = OEField("OEZap_PB_Interaction", Types.FloatVec,
                                     meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapIntEl_field, zapIntEl)

            zapBindSA25_field = OEField("OEZap_SA25_Bind", Types.FloatVec,
                                        meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindSA25_field, zapBindSA25)

            saBuried_field = OEField("OEZap_BuriedArea", Types.FloatVec)
            trajPBSA.set_value(saBuried_field, saBuried)

            zapBindSA6_field = OEField("OEZap_SA6_Bind", Types.FloatVec,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindSA6_field, zapBindSA6)

            # If the OETraj Interaction Energies has been done calculate MMPBSA values
            if 'TrajIntE' in analysesDone:
                opt['Logger'].info('{} found TrajIntE analyses'.format(system_title) )

                # Extract the relevant P-L Interaction Energies from the record
                oeTrjIntERecord = utl.RequestOEFieldType( record, Fields.Analysis.oeintE_rec)
                opt['Logger'].info('{} found TrajIntE record'.format(system_title))

                if self.opt['explicit_water']:

                    PLIntE = utl.RequestOEField(oeTrjIntERecord,
                                                'protein_and_water_ligand_interactionEnergy', Types.FloatVec)
                    opt['Logger'].info('{} found Protein-Water and Ligand force field interaction energies'
                                       .format(system_title))
                else:

                    PLIntE = utl.RequestOEField(oeTrjIntERecord,
                                                'protein_ligand_interactionEnergy', Types.FloatVec)
                    opt['Logger'].info('{} found Protein-Ligand force field interaction energies'
                                       .format(system_title))

                # Calculate  and store MMPB and MMPBSA energies on the trajPBSA record
                zapMMPB = [eInt+eDesol for eInt, eDesol in zip(PLIntE, zapDesolEl)]

                zapMMPB_field = OEField("OEZap_MMPB_Bind", Types.FloatVec,
                                        meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
                trajPBSA.set_value(zapMMPB_field, zapMMPB)

                zapMMPBSA = [eMMPB+eSA6 for eMMPB,eSA6 in zip(zapMMPB, zapBindSA6)]
                zapMMPBSA_field = OEField("OEZap_MMPBSA6_Bind", Types.FloatVec,
                                          meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
                trajPBSA.set_value(Fields.Analysis.zapMMPBSA_fld, zapMMPBSA)

                # Add these energy terms to the earlier dict to store on the record
                PBSAdata['OEZap_MMPB_Bind'] = zapMMPB
                PBSAdata['OEZap_MMPBSA6_Bind'] = zapMMPBSA

                # Clean average MMPBSA to avoid nans and high zap energy values
                avg_mmpbsa, serr_mmpbsa = utl.clean_mean_serr(zapMMPBSA)

                # Add to the record the MMPBSA mean and std
                record.set_value(Fields.Analysis.mmpbsa_traj_mean, avg_mmpbsa)
                record.set_value(Fields.Analysis.mmpbsa_traj_serr, serr_mmpbsa)
                # Add to the record the Average MMPBSA floe report label
                record.set_value(Fields.floe_report_label, "MMPBSA score:<br>{:.1f}  &plusmn; {:.1f} kcal/mol".
                                 format(avg_mmpbsa, serr_mmpbsa))

            # Add the PBSAdata dict to the record
            record.set_value(Fields.Analysis.oepbsa_dict, PBSAdata)

            # Add the trajPBSA record to the parent record
            record.set_value(Fields.Analysis.oepbsa_rec, trajPBSA)
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
        "memory_mb": {"default": 14000},
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
            intE, cplxE, protE, ligE, watE, lwIntE, pwIntE, pw_lIntE = mmpbsa.ProtLigWatInteractionEFromParmedOETraj(
                prmed, ligTraj, protTraj, water_traj, opt)

            if intE is None:
                raise ValueError('{} Calculation of Interaction Energies failed'.format(system_title))

            # protein and ligand traj OEMols now have parmed charges on them; save these
            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ligTraj)
            record.set_value(Fields.Analysis.oetraj_rec, oetrajRecord)

            # make a dict of these energy terms to store on the record
            intEdata = dict()
            intE, cplxE, protE, ligE, watE, lwIntE, pwIntE, pw_lIntE

            # Create new record with traj interaction energy results
            opt['Logger'].info('{} writing trajIntE OERecord'.format(system_title))
            trajIntE = OERecord()

            intE_field = OEField("protein_ligand_interactionEnergy", Types.FloatVec,
                                 meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value(intE_field, intE)

            ligE_field = OEField("ligand_intraEnergy", Types.FloatVec,
                                 meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value(ligE_field, ligE)

            protE_field = OEField("protein_intraEnergy", Types.FloatVec,
                                  meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value(protE_field, protE)

            cplxE_field = OEField("complex_intraEnergy", Types.FloatVec,
                                  meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value(cplxE_field, cplxE)

            watE_field = OEField("water_intraEnergy", Types.FloatVec,
                                 meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))

            trajIntE.set_value(watE_field, watE)

            lwE_field = OEField("ligand_water_interactionEnergy", Types.FloatVec,
                                meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))

            trajIntE.set_value(lwE_field, lwIntE)

            pwE_field = OEField("protein_water_interactionEnergy", Types.FloatVec,
                                meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))

            trajIntE.set_value(pwE_field, pwIntE)

            pw_lIntE_field = OEField("protein_and_water_ligand_interactionEnergy", Types.FloatVec,
                                     meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))

            trajIntE.set_value(pw_lIntE_field, pw_lIntE)

            # Add the intEdata dict to the record
            record.set_value(Fields.Analysis.oeintE_dict, intEdata)

            # Add the trajIntE record to the parent record
            record.set_value(Fields.Analysis.oeintE_rec, trajIntE)
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

    This cube will read in and combine the individual MD traj OEMols for each conformer
    into a single MD traj OEMol for the whole ligand, in preparation for clistering.
    It will do this for both the protein and ligand components of the complex.
    """

    #uuid = "b503c2f4-12e6-49c7-beb6-ee17da177ec2"

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
                .format(system_title) )

            # Go find the ligand and LigTraj fields in each of the conformer records
            if not record.has_field(Fields.Analysis.oetrajconf_rec):
                raise ValueError('{} could not find the conformer record'.format(system_title))
            else:
                opt['Logger'].info('{} found the conformer record'.format(system_title))

            # set up ligand and LigTraj lists then loop over conformer records
            confIdVec = []
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
                confIdVec += [confid]*ligTraj.NumConfs()
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

            if len(ligTrajConfs)<1 or len(protTrajConfs)<1:
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
                system_title, protTraj.NumAtoms(), protTraj.NumConfs()) )


            record.set_value(OEField('ConfIdVec', Types.IntVec), confIdVec)

            # Create new record with OETraj results
            oetrajRecord = OERecord()
            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ligTraj)
            if watTraj:
                oetrajRecord.set_value(OEField('WatTraj', Types.Chem.Mol), watTraj)

            if in_orion():
                collection_id = utl.RequestOEFieldType(record, Fields.collection)
                oetrajRecord.set_value(Fields.collection, collection_id)
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
    This cube gathers together conformers related to the same ligand and their information
    in a new record containing the multi conformer ligand and each conformer record info.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
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

        for sys_id, list_conf_rec in self.lig_sys_ids.items():

            # catch case where for some reason the conf list list_conf_rec is empty
            if len(list_conf_rec)<1:
                print('{} does not have any conformer data'.format(sys_id) )
                continue
            elif len(list_conf_rec)>1:
                # Conformers for each ligand are sorted based on their confid in each ligand record
                list_conf_rec.sort(key=lambda x: x.get_value(Fields.confid))

            new_rec = OERecord()
            new_rec.set_value(Fields.Analysis.oetrajconf_rec, list_conf_rec)
            # Get the first conf to move some general ligand data up to the top level
            rec0 = list_conf_rec[0]
            #   copy all the initial fields in Fields.ligInit_rec up to the top level
            init_rec = rec0.get_value(Fields.ligInit_rec)
            for field in init_rec.get_fields():
                new_rec.set_value(field, init_rec.get_value(field))
            #   next, fields that will simply be copied and not further used here
            protein = rec0.get_value(Fields.protein)
            new_rec.set_value(Fields.protein, protein)
            ligid = rec0.get_value(Fields.ligid)
            new_rec.set_value(Fields.ligid, ligid)
            if in_orion():
                collection_id = rec0.get_value(Fields.collection)
                new_rec.set_value(Fields.collection, collection_id)
            #   finally, fields that will be copied and also further used here
            lig_multi_conf = oechem.OEMol(rec0.get_value(Fields.ligand))
            protein_name = rec0.get_value(Fields.protein_name)

            # if >1 confs, add their confs to the parent ligand at the top level
            for rec in list_conf_rec[1:]:
                lig_multi_conf.NewConf(rec.get_value(Fields.ligand))

            # get name of initial molecule
            init_mol = new_rec.get_value(OEField('Molecule', Types.Chem.Mol))
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


class NMaxWatersLigProt(RecordPortsMixin, ComputeCube):
    title = "NMax Waters"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Ligand', 'Protein', 'Waters']

    description = """
    This cube determines the max number of waters for all the ligands that
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
    # uuid = "a6a11dbb-bc25-4548-bf1a-471bda2f0406"


class ParallelConcatenateTrajMMPBSACube(ParallelMixin, ConcatenateTrajMMPBSACube):
    title = "Parallel " + ConcatenateTrajMMPBSACube.title
    description = "(Parallel) " + ConcatenateTrajMMPBSACube.description
    #uuid = "a6a11dbb-bc25-4548-bf1a-471bda2f0406"


class ParallelTrajPBSACube(ParallelMixin, TrajPBSACube):
    title = "Parallel " + TrajPBSACube.title
    description = "(Parallel) " + TrajPBSACube.description
    uuid = "1f62c9ee-1d83-469c-b6a0-3b20dc64a1e1"


