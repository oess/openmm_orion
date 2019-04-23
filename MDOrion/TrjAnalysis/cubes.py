from cuberecord import OERecordComputeCube

from MDOrion.Standards import Fields

from floereport import FloeReport, LocalFloeReport

from floe.api import parameter

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

from os import environ

import MDOrion.TrjAnalysis.utils as utl

import MDOrion.TrjAnalysis.TrajMMPBSA_utils as mmpbsa

import oetrajanalysis.OETrajBasicAnalysis_utils as oetrjutl

import ensemble2img

from tempfile import TemporaryDirectory


import numpy as np

from openeye import oechem

import oetrajanalysis.Clustering_utils as clusutl

from openeye import oedepict

import os
import traceback

from floe.api import (ParallelMixin)

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
                                                        _clus_floe_report_Trailer,
                                                        trim_svg,
                                                        MakeClusterInfoText)


class MDFloeReportCube(OERecordComputeCube):
    version = "0.1.0"
    title = "MDFloeReportCube"
    description = """
    The floe report cube generates an Orion floe report tiling the input ligands. 
    Each input record must have ligand ID, ligand title, ligand name, the ligand 
    depiction as svg string, the html report string linked to the ligand and 
    optionally the ligand report label. 
    
    Input:
    -------
    Data record Stream - Streamed-in of ligands to be tiled in the Orion floe report

    Output:
    -------
    None
    """
    classification = [["Analysis"]]
    tags = ['Report']

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    upload = parameter.BooleanParameter(
        'upload',
        default=False,
        help_text="Upload floe report to Amazon S3")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.floe_report_dic = dict()

        if in_orion():
            job_id = environ.get('ORION_JOB_ID')
            self.floe_report = FloeReport.start_report("floe_report", job_id=job_id)
        else:
            self.floe_report = LocalFloeReport.start_report("floe_report")

    def process(self, record, port):

        try:

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            system_title = mdrecord.get_title
            system_id = mdrecord.get_id

            if not record.has_value(Fields.floe_report):
                raise ValueError("Missing the report field for the system {}".format(system_title + "_" + system_id))

            report_string = record.get_value(Fields.floe_report)

            if not record.has_value(Fields.ligand_name):
                raise ValueError("Missing the ligand name field")

            ligand_name = record.get_value(Fields.ligand_name)

            if not record.has_value(Fields.floe_report_svg_lig_depiction):
                raise ValueError("Missing the ligand  depiction field")

            ligand_svg = record.get_value(Fields.floe_report_svg_lig_depiction)

            if not record.has_value(Fields.floe_report_label):
                floe_report_label = ""
            else:
                floe_report_label = record.get_value(Fields.floe_report_label)

            self.floe_report_dic[system_id] = (report_string, ligand_svg, ligand_name, floe_report_label)

            # Upload Floe Report
            if self.opt['upload']:

                if in_orion():
                    session = OrionSession()

                    file_upload = File.upload(session,
                                              "{}.html".format(system_title),
                                              report_string)

                    session.tag_resource(file_upload, "floe_report")

                    job_id = environ.get('ORION_JOB_ID')

                    if job_id:
                        session.tag_resource(file_upload, "Job {}".format(job_id))

            self.success.emit(record)

            del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):

        try:
            self.opt['Logger'].info("....Generating Floe Report")

            index = self.floe_report.create_page("index", is_index=True)

            index_content = """
            <style>
            .grid { 
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            grid-gap: 20px;
            align-items: stretch;
            }

            .grid a {
            border: 1px solid #ccc;
            padding: 25px
            }

            .grid svg {
            display: block;  
            max-width: 100%;
            }

            .grid p{
            text-align: center;
            }
            </style>
            <main class="grid">
            """
            # Sort the dictionary keys by using the ligand ID
            for key in sorted(self.floe_report_dic.keys()):

                report_string, ligand_svg, ligand_title, label = self.floe_report_dic[key]

                if len(ligand_title) < 15:
                    page_title = ligand_title
                else:
                    page_title = ligand_title[0:13] + '...'

                page = self.floe_report.create_page(page_title, is_index=False)
                page_link = page.get_link()
                page.set_from_string(report_string)

                index_content += """
                <a href='{}'>
                {}
                <p> {} </p>
                </a>
                """.format(page_link, ligand_svg, label)

            index_content += """
            </main>
            """

            index.set_from_string(index_content)

            self.floe_report.finish_report()

        except Exception as e:
            self.opt['Warning'].warn("It was not possible to generate the floe report: {}".format(str(e)))

        return


class TrajToOEMolCube(ParallelMixin, OERecordComputeCube):
    title = 'Traj to OEMol Cube'

    version = "0.1.0"
    classification = [["Analysis"]]
    tags = ['Trajectory', 'Ligand', 'Protein']

    description = """
    Converting MD Traj into multiconf OEMols for Ligand and Protein
    This cube will take in the MD traj file containing
    the solvated protein:ligand complex and extract
    multiconf OEMols for Ligand and Protein.
    Input parameters:
    -------
    oechem.OEDataRecord - Streamed-in MD results for input
    Output:
    -------
    oechem.OEDataRecord - Stream of output data with trajectory OEMols
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

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
            sys_id = mdrecord.get_id
            opt['Logger'].info('{}: Attempting MD Traj conversion into OEMols'.format(system_title))

            traj_fn = mdrecord.get_stage_trajectory()

            opt['Logger'].info('{} Temp Directory: {}'.format(system_title, os.path.dirname(traj_fn)))
            opt['Logger'].info('{} Trajectory filename: {}'.format(system_title, traj_fn))

            setupOEMol = mdrecord.get_stage_topology(stg_name="System Parametrization")

            opt['Logger'].info('{} Setup topology has {} atoms'.format(system_title, setupOEMol.NumAtoms()))

            # Generate multi-conformer protein and ligand OEMols from the trajectory
            opt['Logger'].info('{} Generating protein and ligand trajectory OEMols'.format(system_title))

            ptraj, ltraj = utl.ExtractAlignedProtLigTraj(setupOEMol, traj_fn)

            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'.format(
                system_title, ptraj.NumAtoms(), ptraj.NumConfs()))
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'.format(
                system_title, ltraj.NumAtoms(), ltraj.NumConfs()))

            # Generate average and median protein and ligand OEMols from ptraj, ltraj
            opt['Logger'].info('{} Generating protein and ligand median and average OEMols'.format(system_title))

            ligMedian, protMedian, ligAverage, protAverage = oetrjutl.AnalyseProteinLigandTrajectoryOEMols(ltraj, ptraj)

            # Generate interactive trajectory SVG
            opt['Logger'].info('{} Generating interactive trajectory SVG'.format(system_title))
            trajSVG = ensemble2img.run_ensemble2img(ligMedian, protMedian, ltraj, ptraj)


            # Overwrite MDStages with only first (setup) and last (production) stages
            # newMDStages = [md_stage0_record, md_stageLast_record]
            # record.set_value(Fields.md_stages, newMDStages)


            # Create new record with OETraj results
            oetrajRecord = OERecord()

            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ltraj)

            oetrajRecord.set_value(OEField('LigMedian', Types.Chem.Mol), ligMedian)

            oetrajRecord.set_value(OEField('ProtMedian', Types.Chem.Mol), protMedian)

            oetrajRecord.set_value(OEField('LigAverage', Types.Chem.Mol), ligAverage)

            oetrajRecord.set_value(OEField('ProtAverage', Types.Chem.Mol), protAverage)

            if in_orion():
                oetrajRecord.set_value(Fields.collection, mdrecord.collection_id)

            mdrecord_traj = MDDataRecord(oetrajRecord)

            mdrecord_traj.set_protein_traj(ptraj, shard_name="ProteinTrajConfs_" + system_title + '_' + str(sys_id))

            TrajSVG_field = OEField('TrajSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))

            oetrajRecord.set_value(TrajSVG_field, trajSVG)

            record.set_value(OEField('OETraj', Types.Record), oetrajRecord)

            # update or initiate the list of analyses that have been done
            analysesDoneField = OEField('AnalysesDone', Types.StringVec)

            if record.has_value(analysesDoneField):
                analysesDone = utl.RequestOEFieldType(record, analysesDoneField)
                analysesDone.append('OETraj')
            else:
                analysesDone = ['OETraj']

            record.set_value( analysesDoneField, analysesDone)
            opt['Logger'].info('{}: saved protein and ligand traj OEMols'.format(system_title))

            self.success.emit(record)

            del mdrecord
            del mdrecord_traj

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class TrajPBSACube(ParallelMixin, OERecordComputeCube):
    title = "Trajectory Poisson-Boltzmann and Surface Area Energies"
    version = "0.0.0"
    classification = [["Analysis"]]
    tags = ['OEChem', 'Zap', 'TrajAnalysis', 'MMPBSA']
    description = """
    Protein-ligand interaction solvation energies are calculated on an existing MD trajectory.
    The trajectory is taken from pre-existing protein and ligand trajectory OEMols.
    The Poisson-Boltzmann and Surface Area methods in the OEZap toolkits are  used.

    Input:
    -------
    Data Record with the ligand and Protein trajectory OEMols; the trajectory frames are
    included as conformers on the molecule.

    Output:
    -------
    Data Record - The various energy components associated with the Poisson-Boltzmann and
    Surface Area energies are attached to the record as per-frame vectors of floats.
    The energy units are in kcal/mol.
    """

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
            opt['Logger'].info(' Beginning TrajPBSACube')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to compute MD Traj PBSA energies'.format(system_title))

            # Check that the OETraj analysis has been done
            analysesDone = utl.RequestOEField( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title))
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title))

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = utl.RequestOEField(record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title))
            ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                               .format(system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()))

            mdtrajrecord = MDDataRecord(oetrajRecord)
            protTraj = mdtrajrecord.get_protein_traj

            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                               .format(system_title, protTraj.NumAtoms(), protTraj.NumConfs()))

            # Compute PBSA energies for the protein-ligand complex
            zapBind, zapBindPB, zapDesolEl, zapIntEl, zapBindSA25, saBuried = mmpbsa.TrajPBSA(
                                       ligTraj, protTraj)
            if zapBind is None:
                raise ValueError('{} Calculation of PBSA energies failed'.format(system_title) )
            # generate Surface Areas energy for buried SA based on 0.006 kcal/mol/A^2
            zapBindSA6 = [sa*-0.006 for sa in saBuried]

            # Create new record with traj interaction energy results
            opt['Logger'].info('{} writing trajPBSA OERecord'.format(system_title) )
            trajPBSA = OERecord()
            #
            zapBind_field = OEField("OEZap_PBSA25_Bind", Types.FloatVec,
                                    meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBind_field, zapBind)
            #
            zapBindPB_field = OEField("OEZap_PB_Bind", Types.FloatVec,
                                      meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindPB_field, zapBindPB)
            #
            zapDesolEl_field = OEField("OEZap_PB_Desolvation", Types.FloatVec,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapDesolEl_field, zapDesolEl)
            #
            zapIntEl_field = OEField("OEZap_PB_Interaction", Types.FloatVec,
                                     meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapIntEl_field, zapIntEl)
            #
            zapBindSA25_field = OEField("OEZap_SA25_Bind", Types.FloatVec,
                                        meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindSA25_field, zapBindSA25)
            #
            saBuried_field = OEField("OEZap_BuriedArea", Types.FloatVec)
            trajPBSA.set_value(saBuried_field, saBuried)
            #
            zapBindSA6_field = OEField("OEZap_SA6_Bind", Types.FloatVec,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindSA6_field, zapBindSA6)

            # If the OETraj Interaction Energies has been done calculate MMPBSA values
            if 'TrajIntE' in analysesDone:
                opt['Logger'].info('{} found TrajIntE analyses'.format(system_title) )

                # Extract the relevant P-L Interaction Energies from the record
                oeTrjIntERecord = utl.RequestOEField(record, 'TrajIntE', Types.Record)
                opt['Logger'].info('{} found TrajIntE record'.format(system_title))
                PLIntE = utl.RequestOEField(oeTrjIntERecord,
                                            'protein_ligand_interactionEnergy', Types.FloatVec)
                opt['Logger'].info('{} found Protein-Ligand force field interaction energies'
                                   .format(system_title))
                # Calculate  and store MMPB and MMPBSA energies on the trajPBSA record
                zapMMPB = [eInt+eDesol for eInt,eDesol in zip(PLIntE, zapDesolEl)]
                zapMMPB_field = OEField("OEZap_MMPB_Bind", Types.FloatVec,
                                        meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
                trajPBSA.set_value(zapMMPB_field, zapMMPB)
                #
                zapMMPBSA = [eMMPB+eSA6 for eMMPB,eSA6 in zip(zapMMPB, zapBindSA6)]
                zapMMPBSA_field = OEField("OEZap_MMPBSA6_Bind", Types.FloatVec,
                                          meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
                trajPBSA.set_value(zapMMPBSA_field, zapMMPBSA)

                # Average MMPBSA for all the trajectory frames
                np_mmpbsa = np.array(zapMMPBSA)
                avg_mmpbsa = np_mmpbsa.mean()
                std_mmpbsa = np_mmpbsa.std()

                # Add to the record the Average MMPBSA floe report label
                record.set_value(Fields.floe_report_label, "MMPBSA = {:.1f}  &plusmn; {:.1f} kcal/mol".
                                 format(avg_mmpbsa, std_mmpbsa))

            # Add the trajPBSA record to the parent record
            record.set_value(OEField('TrajPBSA', Types.Record), trajPBSA)
            analysesDone.append('TrajPBSA')
            record.set_value(OEField('AnalysesDone', Types.StringVec), analysesDone)
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


class TrajInteractionEnergyCube(ParallelMixin, OERecordComputeCube):
    title = "Trajectory Interaction Energies"
    version = "0.0.0"
    classification = [["Analysis"]]
    tags = ['OEChem', 'OpenMM', 'TrajAnalysis', 'MMPBSA']
    description = """
    Protein-ligand interaction energies are calculated on an existing MD trajectory.
    The trajectory is taken from pre-existing protein and ligand trajectory OEMols.
    The forcefield used is taken from the parmed object associated with the trajectory
    OEMols.

    Input:
    -------
    Data Record with the ligand and Protein trajectory OEMols; the trajectory frames are
    included as conformers on the molecule. The associated parmed object is also on the
    record.

    Output:
    -------
    Data Record - The various energy components associated with the protein-ligand
    interaction energies are attached to the record as per-frame vectors of floats.
    This includes the MM interaction potential energies and their ligand, protein
    and complex components. The energy units are in kcal/mol.
    """

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
            opt['Logger'].info(' Beginning TrajInteractionEnergyCube')

            mdrecord = MDDataRecord(record)

            system_title = mdrecord.get_title

            opt['Logger'].info('{} Attempting to compute MD Traj protein-ligand Interaction energies'
                .format(system_title) )

            # Check that the OETraj analysis has been done
            analysesDone = utl.RequestOEField( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = utl.RequestOEField( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                .format( system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )

            mdtrajrecord = MDDataRecord(oetrajRecord)
            protTraj = mdtrajrecord.get_protein_traj

            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                .format( system_title, protTraj.NumAtoms(), protTraj.NumConfs()) )

            prmed = mdrecord.get_parmed(sync_stage_name='last')

            # Compute interaction energies for the protein, ligand, and complex subsystems
            intE, cplxE, protE, ligE = mmpbsa.ProtLigInteractionEFromParmedOETraj(
                                       prmed, ligTraj, protTraj)
            if intE is None:
                raise ValueError('{} Calculation of Interaction Energies failed'.format(system_title) )

            # protein and ligand traj OEMols now have parmed charges on them; save these
            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ligTraj)

            record.set_value(OEField('OETraj', Types.Record), oetrajRecord)

            # Create new record with traj interaction energy results
            opt['Logger'].info('{} writing trajIntE OERecord'.format(system_title))
            trajIntE = OERecord()
            #
            intE_field = OEField("protein_ligand_interactionEnergy", Types.FloatVec,
                                meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( intE_field, intE)
            #
            ligE_field = OEField("ligand_intraEnergy", Types.FloatVec,
                                meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( ligE_field, ligE)
            #
            protE_field = OEField("protein_intraEnergy", Types.FloatVec,
                                 meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( protE_field, protE)
            #
            cplxE_field = OEField("complex_intraEnergy", Types.FloatVec,
                                 meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( cplxE_field, cplxE)
            # Add the trajIntE record to the parent record
            record.set_value( OEField( 'TrajIntE', Types.Record), trajIntE)
            analysesDone.append( 'TrajIntE')
            record.set_value( OEField( 'AnalysesDone', Types.StringVec), analysesDone)
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


class ClusterOETrajCube(ParallelMixin, OERecordComputeCube):
    title = 'Cluster Protein-Ligand Traj OEMols'

    version = "0.1.0"
    classification = [["Analysis"]]
    tags = ['Clustering', 'Ligand', 'Protein']

    description = """
    Cluster  multiconf MD trajectory OEMols for Ligand and Protein

    This cube will take in the MD traj OEMols containing
    the protein and ligand components of the complex and cluster
    them based on ligand RMSD.

    Input parameters:
    -------
    oechem.OEDataRecord - Streamed-in input data for the system to cluster

    Output:
    -------
    oechem.OEDataRecord - Stream of output data for the clustered system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
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
            opt['Logger'].info(' Beginning ClusterOETrajCube')
            system_title = utl.RequestOEFieldType( record, Fields.title)
            opt['Logger'].info('{} Attempting to cluster MD Traj'
                .format(system_title) )

            # Check that the OETraj analysis has been done
            analysesDone = utl.RequestOEField( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = utl.RequestOEField( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                .format( system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )

            mdtrajrecord = MDDataRecord(oetrajRecord)
            protTraj = mdtrajrecord.get_protein_traj

            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                .format( system_title, protTraj.NumAtoms(), protTraj.NumConfs()) )

            # Cluster ligand traj into cluster OEMols with matching protein OEMols
            opt['Logger'].info('{} starting clustering'.format(system_title) )
            clusResults = clusutl.ClusterLigTraj( ligTraj)
            opt['Logger'].info('{} plotting RMSD histogram'.format(system_title) )
            trajHistRMSD_svg = clusutl.ClusterLigTrajHistRMSD(clusResults)
            # Heat map is large, time-consuming, and non-essential... put back if needed
            #opt['Logger'].info('{} plotting RMSD heat map'.format(system_title) )
            #trajHeatRMSD_png = clusutl.ClusterLigTrajHeatRMSD(clusResults)
            opt['Logger'].info('{} plotting cluster strip plot'.format(system_title) )
            trajClus_svg = clusutl.ClusterLigTrajClusPlot(clusResults)
            # Calculate RMSD of ligand traj from ligand initial pose
            vecRmsd = oechem.OEDoubleArray( ligTraj.GetMaxConfIdx() )
            ligInitPose = utl.RequestOEFieldType( record, Fields.ligand)
            oechem.OERMSD( ligInitPose, ligTraj, vecRmsd)
            rmsdInit_svg = clusutl.RmsdFromInitialPosePlot( clusResults['ClusterVec'], vecRmsd)

            # Create trajSVG for each major cluster (major= 10% or more of traj)
            clusLigAvgMol = []
            clusProtAvgMol = []
            clusTrajSVG = []
            for clusID, count in enumerate(clusResults['ClusterCounts']):
                if clusResults['nFrames']/count<=10:
                    opt['Logger'].info('generating svg for {} cluster {}'.format( system_title, clusID))
                    clusLig = clusutl.TrajOEMolFromCluster( ligTraj, clusResults['ClusterVec'], clusID)
                    opt['Logger'].info( 'new mol {} with {} confs:'.format( clusLig.GetTitle(),clusLig.NumConfs()) )
                    clusProt = clusutl.TrajOEMolFromCluster( protTraj, clusResults['ClusterVec'], clusID)
                    opt['Logger'].info( 'new mol {} with {} confs:'.format( clusProt.GetTitle(),clusProt.NumConfs()) )
                    ligMed, protMed, ligAvg, protAvg = oetrjutl.AnalyseProteinLigandTrajectoryOEMols( clusLig, clusProt)
                    clusSVG = ensemble2img.run_ensemble2img(ligAvg, protAvg, clusLig, clusProt)
                    clusLigAvgMol.append( ligAvg)
                    clusProtAvgMol.append( protAvg)
                    clusTrajSVG.append( clusSVG)

            # Create new record with trajClus results
            opt['Logger'].info('{} writing trajClus OERecord'.format(system_title) )
            trajClus = OERecord()
            #
            ClusLigAvgMol_field = OEField( 'ClusLigAvgMol', Types.Chem.MolVec)
            trajClus.set_value( ClusLigAvgMol_field, clusLigAvgMol)
            #
            ClusProtAvgMol_field = OEField( 'ClusProtAvgMol', Types.Chem.MolVec)
            trajClus.set_value( ClusProtAvgMol_field, clusProtAvgMol)
            #
            ClusTrajSVG_field = OEField( 'ClusTrajSVG', Types.StringVec)
            trajClus.set_value( ClusTrajSVG_field, clusTrajSVG)
            #
            ClusterMethod_field = OEField( 'ClusterMethod', Types.String)
            trajClus.set_value( ClusterMethod_field, clusResults['ClusterMethod'])
            #
            HDBSCAN_alpha_field = OEField( 'HDBSCAN_alpha', Types.Float)
            trajClus.set_value( HDBSCAN_alpha_field, clusResults['HDBSCAN_alpha'])
            #
            nFrames_field = OEField( 'nFrames', Types.Int)
            trajClus.set_value( nFrames_field, clusResults['nFrames'])
            #
            nClusters_field = OEField( 'nClusters', Types.Int)
            trajClus.set_value( nClusters_field, clusResults['nClusters'])
            #
            clusterCounts_field = OEField( 'ClusterCounts', Types.IntVec)
            trajClus.set_value( clusterCounts_field, clusResults['ClusterCounts'])
            #
            Clusters_field = OEField( 'Clusters', Types.IntVec)
            trajClus.set_value( Clusters_field, clusResults['ClusterVec'])
            #
            HistSVG_field = OEField( 'HistSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            trajClus.set_value( HistSVG_field, trajHistRMSD_svg)
            #
            ClusSVG_field = OEField( 'ClusSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            trajClus.set_value( ClusSVG_field, trajClus_svg)
            # Heat map is large, time-consuming, and non-essential... put back if needed
            # HeatPNG_field = OEField( 'HeatPNG', Types.Blob, meta=OEFieldMeta().set_option(Meta.Hints.Image_PNG))
            # trajClus.set_value( HeatPNG_field, trajHeatRMSD_png)
            #
            rmsdInit_field = OEField( 'rmsdInitPose', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            trajClus.set_value( rmsdInit_field, rmsdInit_svg)

            #
            record.set_value( OEField( 'TrajClus', Types.Record), trajClus)
            analysesDone.append( 'TrajClus')
            record.set_value( OEField( 'AnalysesDone', Types.StringVec), analysesDone)
            opt['Logger'].info('{} finished writing trajClus OERecord'.format(system_title) )

            self.success.emit(record)

            del mdtrajrecord

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ClusterOETrajCube on {}'.format(str(e),system_title) )
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class MDTrajAnalysisClusterReport(ParallelMixin, OERecordComputeCube):
    title = 'Extract relevant outputs of MD Traj Cluster  Analysis'

    version = "0.1.0"
    classification = [["Analysis"]]
    tags = ['Ligand', 'Protein']

    description = """
    Extract relevant outputs of Ligand and Protein
    Short Traj MD Traj Analysis and write them to files.

    This cube takes as input the OERecord containing the work
    product of trajectory analysis on Short Traj MD results.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
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

            # title of entire solvated protein-ligand system
            opt['Logger'].info('Starting Floe Report generation for MD Traj Analysis')

            system_title = utl.RequestOEFieldType(record, Fields.title)

            opt['Logger'].info('{} Attempting to extract MD Traj Analysis results'.format(system_title))

            ligInitPose = utl.RequestOEFieldType(record, Fields.ligand)

            lig_name = utl.RequestOEFieldType(record, Fields.ligand_name)

            protInitPose = utl.RequestOEFieldType(record, Fields.protein)

            asiteSVG = utl.PoseInteractionsSVG(ligInitPose, protInitPose, width=400, height=265)

            # Extract the traj SVG from the OETraj record
            analysesDone = utl.RequestOEField(record, 'AnalysesDone', Types.StringVec)

            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )

            # Extract the relevant traj SVG from the OETraj record
            oetrajRecord = utl.RequestOEField(record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title))

            trajSVG = utl.RequestOEField(oetrajRecord, 'TrajSVG', Types.String)

            # Extract the label for the MMPBSA score for the whole trajectory
            mmpbsaLabelStr = utl.RequestOEField(record, 'Floe_report_label_OPLMD', Types.String)

            # Extract Ligand average Bfactor
            ligand_bfactor = utl.RequestOEField(oetrajRecord, 'LigAverage', Types.Chem.Mol)

            # Extract the three plots from the TrajClus record
            analysesDone = utl.RequestOEField(record, 'AnalysesDone', Types.StringVec)

            if 'TrajClus' not in analysesDone:
                raise ValueError('{} does not have TrajClus analyses done'.format(system_title))
            else:
                opt['Logger'].info('{} found TrajClus analyses'.format(system_title))

            # Extract the relevant traj SVG from the TrajClus record
            clusRecord = utl.RequestOEField(record, 'TrajClus', Types.Record)

            opt['Logger'].info('{} found TrajClus record'.format(system_title) )
            trajHistRMSD_svg = utl.RequestOEField(clusRecord, 'HistSVG', Types.String)
            trajClus_svg = utl.RequestOEField(clusRecord, 'ClusSVG', Types.String)
            rmsdInit_svg = utl.RequestOEField(clusRecord, 'rmsdInitPose', Types.String)
            clusTrajSVG = utl.RequestOEField(clusRecord, 'ClusTrajSVG', Types.StringVec)

            opt['Logger'].info('{} found the TrajClus plots'.format(system_title))

            # Generate text string about Clustering information
            clusData = {}

            clusData['nFrames'] = utl.RequestOEField(clusRecord, 'nFrames', Types.Int)
            clusData['ClusterMethod'] = utl.RequestOEField(clusRecord, 'ClusterMethod', Types.String)
            clusData['HDBSCAN_alpha'] = utl.RequestOEField(clusRecord, 'HDBSCAN_alpha', Types.Float)
            clusData['nClusters'] = utl.RequestOEField(clusRecord, 'nClusters', Types.Int)
            clusData['ClusterVec'] = utl.RequestOEField(clusRecord, 'Clusters', Types.IntVec)
            clusData['ClusterCounts'] = utl.RequestOEField(clusRecord, 'ClusterCounts', Types.IntVec)

            opt['Logger'].info('{} finished writing analysis files'.format(system_title))

            # prepare the 2D structure depiction
            oedepict.OEPrepareDepiction(ligInitPose)
            img = oedepict.OEImage(400, 300)
            oedepict.OERenderMolecule(img, ligInitPose)

            # get the palette of graph marker colors
            nClustersP1 = clusData['nClusters']+1
            clusRGB = utl.ColorblindRGBMarkerColors(nClustersP1)
            clusRGB[-1] = (76, 76, 76)

            with TemporaryDirectory() as output_directory:

                # write the report
                reportFName = os.path.join(output_directory, system_title + '_ClusReport.html')

                report_file = open(reportFName, 'w')

                report_file.write(_clus_floe_report_header)

                for i in range(len(clusTrajSVG)+2):
                    report_file.write("""
                  div.cb-floe-report__tab-wrapper input:nth-of-type({clusID}):checked ~ .cb-floe-report__tab-content:nth-of-type({clusID}) {{ display: block; }}
                """.format(clusID=i+1))

                report_file.write(_clus_floe_report_header2)

                report_file.write(_clus_floe_report_midHtml0.format(
                    query_depiction=oedepict.OEWriteImageToString("svg", img).decode("utf-8")))

                report_file.write("""      <h3>
                        {mmpbsaLabel}
                      </h3>""".format(mmpbsaLabel=mmpbsaLabelStr ))

                analysis_txt = MakeClusterInfoText(clusData, clusRGB)
                report_file.write("".join(analysis_txt))

                report_file.write(_clus_floe_report_midHtml1)

                report_file.write("""      <input type="radio" name="tab" id="cb-floe-report__tab-1-header" checked>
                      <label class="cb-floe-report__tab-label" for="cb-floe-report__tab-1-header">Overall</label>""")

                CurrentTabId = 1

                for i, (clus, rgb) in enumerate(zip(clusTrajSVG, clusRGB)):
                    CurrentTabId = i+2
                    report_file.write("""      <input type="radio" name="tab" id="cb-floe-report__tab-{tabID}-header">
                      <label class="cb-floe-report__tab-label" for="cb-floe-report__tab-{tabID}-header" style="
                                background-color: rgb({r},{g},{b});
                                color: white;">Cluster {clusNum}</label>
                                """.format(tabID=CurrentTabId, clusNum=i, r=rgb[0], g=rgb[1], b=rgb[2]))

                report_file.write("""      <input type="radio" name="tab" id="cb-floe-report__tab-{tabID}-header">
                      <label class="cb-floe-report__tab-label" for="cb-floe-report__tab-{tabID}-header">Initial Pose</label>
                      """.format(tabID=CurrentTabId+1, clusNum=i ))

                report_file.write("""      <div class="cb-floe-report__tab-content">
                        {traj}
                      </div>""".format(traj=trim_svg(trajSVG)))

                for clusSVG in clusTrajSVG:
                    report_file.write("""      <div class="cb-floe-report__tab-content">
                        {traj}
                      </div>
                      """.format(traj=trim_svg(clusSVG)))

                report_file.write("""      <div class="cb-floe-report__tab-content">
                        {traj}
                      </div>
                      """.format(traj=trim_svg(asiteSVG)))

                report_file.write(_clus_floe_report_midHtml2)

                report_file.write(_clus_floe_report_Trailer.format(
                    clusters=trim_svg(trajClus_svg),
                    rmsdInit=trim_svg(rmsdInit_svg)))

                report_file.close()

                with open(reportFName, 'r') as f:
                    report_html_str = f.read()

                record.set_value(Fields.floe_report, report_html_str)

                # Copy Bfactors from the average Bfactor ligand to a copy of the ligand initial pose
                ligand_init = oechem.OEMol(ligInitPose)

                for at_avg_bfac, at_init in zip(ligand_bfactor.GetAtoms(), ligand_init.GetAtoms()):
                    if at_avg_bfac.GetAtomicNum() == at_init.GetAtomicNum():
                        res_avg_bfac = oechem.OEAtomGetResidue(at_avg_bfac)
                        bfactor_avg = res_avg_bfac.GetBFactor()
                        res_init = oechem.OEAtomGetResidue(at_init)
                        res_init.SetBFactor(bfactor_avg)
                        oechem.OEAtomSetResidue(at_init, res_init)
                    else:
                        raise ValueError("Atomic number mismatch {} vs {}".format(at_avg_bfac.GetAtomicNum(),
                                                                                  at_init.GetAtomicNum()))
                # Create svg for the report tile
                lig_svg = utl.ligand_to_svg_stmd(ligand_init, lig_name)

                record.set_value(Fields.floe_report_svg_lig_depiction, lig_svg)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


# import traceback
#
# from floe.api import ParallelMixin, parameter
#
# from cuberecord import OERecordComputeCube
#
# from cuberecord.ports import RecordInputPort
#
# from datarecord import (OEField,
#                         Types,
#                         OERecord)
#
#
# from tempfile import TemporaryDirectory
#
# import os
#
# from oeommtools.utils import split
#
# from Standards import Fields, MDStageTypes
#
# import mdtraj as md
#
# import sstmap as sm
#
# from MDEngines import utils as omm_utils
#
# from floe.constants import *
#
# from openeye import oechem, oegrid
#
# from TrjAnalysis import sstmap_utils
#
# from TrjAnalysis.sstmap_utils import GISTFields
#
# import copy as cp
#
# import shutil
#
# from orionclient.session import in_orion
#
# from shutil import copyfile
#
#
# class SSTMapHsa(ParallelMixin, OERecordComputeCube):
#
#     version = "0.1.0"
#
#     title = "SSTMAP HSA Analysis"
#
#     description = """
#         SSTMap performs Water Thermodynamics analysis.
#         SSTMaps supports hydration site analysis (HSA)
#         and Grid Inhomogeneous Solvation Theory (GIST).
#
#         SSTMap has been developed at Kurtzman Lab Lehman College
#         For more details, please visit
#         sstmap.org @ https://github.com/KurtzmanLab/SSTMap
#         """
#     classifications = [["SSTMap Analysis", "SSTMap HSA"]]
#
#     tags = [tag for lists in classifications for tag in lists]
#
#     # Override defaults for some parameters
#     parameter_overrides = {
#         "memory_mb": {"default": 6000},
#         "spot_policy": {"default": "Allowed"},
#         "prefetch_count": {"default": 1},  # 1 molecule at a time
#         "item_count": {"default": 1}  # 1 molecule at a time
#     }
#
#     start_frame = parameter.IntegerParameter(
#         'start_frame',
#         default=0,
#         min_value=0,
#         level=ADVANCED,
#         help_text="Frame index to start the SSTMap analysis. Default: 0."
#     )
#
#     total_frames = parameter.IntegerParameter(
#         'total_frames',
#         max_value=100000,
#         default=100,
#         level=ADVANCED,
#         help_text="Total number of frames to process during the analysis. Default: 100."
#     )
#
#     hsa_rad = parameter.DecimalParameter(
#         'hsa_rad',
#         min_value=5.0,
#         max_value=10.0,
#         default=5.0,
#         help_text="Distance cutoff (in Angstrom) used to identify hsa region. All waters within this distance from any"
#                   " of the ligand atom are included in the analysis. Default: 5.0."
#     )
#
#     lig_res_name = parameter.StringParameter(
#         'lig_res_name',
#         default='LIG',
#         max_length=4,
#         help_text="Resname to use to identify the ligand"
#     )
#
#     # Values taken from AMBER 17 manual p. 610
#     # Water Model    Mean Energy (Eww-norm) (kcal/mol/water)    Number Density (A^-3)
#     # TIP3P          -9.533                                     0.0329
#     # TIP4PEW        -11.036                                    0.0332
#     # TIP4P          -9.856                                     0.0332
#     # TIP5P          -9.596                                     0.0329
#     # TIP3PFW        -11.369                                    0.0334
#     # SPCE           -11.123                                    0.0333
#     # SPCFW          -11.873                                    0.0329
#     # OPC                                                       0.0333
#
#     wat_model=parameter.StringParameter(
#         'wat_model',
#         default='TIP3P',
#         choices=['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC'],
#         level=ADVANCED,
#         help_text="Water model used during the simulation. Used to set bulk density number. Default: TIP3P."
#     )
#
#     ligand_port = RecordInputPort("ligand_port", initializer=True)
#
#     # Uncomment this and implement if you need to initialize the cube
#     def begin(self):
#         self.opt = vars(self.args)
#         self.opt['Logger'] = self.log
#
#         # Generate dictionary of water models and bulk density
#         wat_model = ['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC']
#         wat_model_bulk_density = [0.0329, 0.0332, 0.0332, 0.0329, 0.0334, 0.0333, 0.0329, 0.0333]
#         self.wat_model_density_dic = dict(zip(wat_model,wat_model_bulk_density))
#
#         for ligand in self.ligand_port:
#             self.ligand = ligand.get_value(Fields.primary_molecule)
#
#     # Records are passed to this function for processing.
#     def process(self, record, port):
#         try:
#             if port == "intake":
#
#                 # Generation options dictionary
#                 opt = dict(self.opt)
#
#                 # Assigning the bulk density bas on water model
#                 opt['rho_bulk'] = self.wat_model_density_dic[opt['wat_model']]
#
#                 if not record.has_value(Fields.primary_molecule):
#                     self.log.error("Missing molecule Primary Molecule' field")
#                     self.failure.emit(record)
#                     return
#
#                 system = record.get_value(Fields.primary_molecule)
#
#                 if not record.has_value(Fields.title):
#                     self.log.warn("Missing record Title field")
#                     system_title = system.GetTitle()[0:12]
#                 else:
#                     system_title = record.get_value(Fields.title)
#
#                 if not record.has_value(Fields.id):
#                     raise ValueError("Missing ID Field")
#
#                 sys_id = record.get_value(Fields.id)
#
#                 sys_info = system_title + '_' + str(sys_id)
#
#                 # Get the MDStageRecord list from the record
#                 if record.has_value(Fields.md_stages):
#                     mdstages = record.get_value(Fields.md_stages)
#                 else:
#                     raise ValueError("Field md_stages is missing!")
#
#                 # Get the MDStageRecord for the production stage from the MDStageRecord list
#                 # That correspond to the last member of the MDStageRecord list
#                 mdstage_prod = mdstages[-1]
#
#                 # Get the MDSystemRecord for the production stage
#                 if mdstage_prod.has_value(Fields.md_system):
#                     mdsystem_prod = mdstage_prod.get_value(Fields.md_system)
#                 else:
#                     raise ValueError("Field md_system is missing!")
#
#                 # Get the PARMED object from the MDSystemRecord
#                 prod_topology_parmed = mdsystem_prod.get_value(Fields.structure)
#
#                 # Generate the parameter supplementary file using the parmed object
#                 with TemporaryDirectory() as output_directory:
#                     opt['Logger'].info("{} - Temporary directory: {}\n".format(sys_info, output_directory))
#
#                     # Get the name of the trajectory from the  production MDStageRecord
#                     if mdstage_prod.has_value(Fields.trajectory):
#
#                         if in_orion():
#                             prod_traj_path = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                             prod_traj_filename = os.path.join(output_directory, "trajectory.h5")
#                             copyfile(prod_traj_path, prod_traj_filename)
#                         else:
#                             prod_traj_filename = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                     else:
#                         raise ValueError("MD_stages do not have a trajectory!")
#
#                     opt['Logger'].info("{} - Trajectory file name: {}".format(sys_info, prod_traj_filename))
#
#                     # Get the final structure from the production stage
#                     prod_coord_eomol = mdsystem_prod.get_value(Fields.topology)
#
#                     # Extract the ligand from the final frame
#                     prot, lig, wat, excp = split(prod_coord_eomol)
#
#                     self.log.info("System name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
#                                   "Water atom numbers = {}\nExcipients atom numbers = {}".format(sys_info,
#                                                                                                  prot.NumAtoms(),
#                                                                                                  lig.NumAtoms(),
#                                                                                                  wat.NumAtoms(),
#                                                                                                  excp.NumAtoms()))
#
#                     # Generate Variables needed for running SSTMap
#                     ligand_filename = os.path.join(output_directory, "ligand.pdb")
#
#                     ofs = oechem.oemolostream(ligand_filename)
#                     pdb_flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
#                     ofs.SetFlavor(oechem.OEFormat_PDB, pdb_flavor)
#                     if lig.GetMaxAtomIdx() > 0:
#                         oechem.OEWriteConstMolecule(ofs, lig)
#                         ligand_align = False
#                     else:
#                         oechem.OEWriteConstMolecule(ofs, self.ligand)
#                         ligand_align = True
#                     ofs.close()
#
#                     if ligand_align:
#                         mdstage_setup = mdstages[0]
#                         if mdstage_setup.get_value(Fields.stage_type) == MDStageTypes.SETUP:
#                             mdsytem_setup = mdstage_setup.get_value(Fields.md_system)
#                             setup_topology = mdsytem_setup.get_value(Fields.topology)
#
#                     opt['Logger'].info("{} - Processing Trajectory".format(sys_info))
#                     # In order to get meaningful energy values we need to strip the ions from the trajectory
#                     # Also, we need to fit the protein to a reference frame to remove translations and
#                     # rotations.
#                     top_filename, parm_filename, prod_traj_filename, aligned_prot_oemol = sstmap_utils.\
#                         process_trajectory(prod_traj_filename,
#                                            prod_topology_parmed,
#                                            opt['lig_res_name'],
#                                            output_directory,
#                                            reference_topology=setup_topology)
#
#                     # Get number of frames
#                     total_number_frames = 0
#                     if prod_traj_filename.endswith(".h5"):
#                         for chunk in md.iterload(prod_traj_filename):
#                             total_number_frames += chunk.n_frames
#                     else:
#                         for chunk in md.iterload(prod_traj_filename, top=top_filename):
#                             total_number_frames += chunk.n_frames
#
#                     # Start SSTMap HSA analysis
#                     # Change to tmp directory to avoid data overwrite
#                     cwd = os.getcwd()
#                     os.chdir(output_directory)
#
#                     opt['Logger'].info("{} - Starting HSA Calculation....".format(sys_info))
#
#                     # Initialize HSA calculation
#                     hsa = sm.SiteWaterAnalysis(topology_file=top_filename,
#                                                trajectory=prod_traj_filename,
#                                                start_frame=0,
#                                                num_frames=total_number_frames,
#                                                supporting_file=parm_filename,
#                                                ligand_file=ligand_filename,
#                                                hsa_region_radius=opt['hsa_rad'])
#
#                     # Initialize hydration sites
#                     hsa.initialize_hydration_sites()
#
#                     # Print System summary
#                     hsa.print_system_summary()
#
#                     # Get frame information
#                     cluster_frame_info_list = cp.deepcopy(hsa.site_waters)
#
#                     # Generate clusters and calculate quantities
#                     hsa.calculate_site_quantities()
#
#                     # Write Calculation summary
#                     hsa.write_calculation_summary()
#
#                     # Write data
#                     hsa.write_data()
#                     ############
#
#                     # Generate EOMol for Water cluster depiction
#                     multi_confomer_cluster_list = sstmap_utils.process_clusters(output_directory, total_number_frames, cluster_frame_info_list)
#
#                     # Generate OEMOL for Most probable configuration
#                     most_prob_config = sstmap_utils.probable_conf(output_directory)
#
#                     # Create new record with results
#                     new_record = OERecord()
#                     new_record.set_value(Fields.primary_molecule, most_prob_config)
#
#                     new_record.set_value(Fields.protein, aligned_prot_oemol)
#
#                     # Create Field
#                     mol_vec = OEField("clusters", Types.Chem.MolVec)
#                     new_record.set_value(mol_vec, multi_confomer_cluster_list)
#
#                     alloutfile = os.path.join(output_directory, 'MC_all_clusters.oeb')
#                     allofs = oechem.oemolostream(alloutfile)
#
#                     for conf in multi_confomer_cluster_list:
#                         oechem.OEWriteConstMolecule(allofs, conf)
#
#                     allofs.close()
#
#                     hsa_data = os.path.join(cwd, "HSA_Results_data")
#                     shutil.copytree(output_directory, hsa_data)
#
#                     self.success.emit(new_record)
#
#         except:
#             # Attach an error message to the molecule that failed
#             self.log.error(traceback.format_exc())
#             # Return failed mol
#             self.failure.emit(record)
#
#
# class SSTMapGist(ParallelMixin, OERecordComputeCube):
#     # Cube documentation.  This documentation for this cube, and all other cubes in this repository, can be converted
#     # to html by calling 'invoke docs' from the root directory of this repository.  This documentation will also
#     # appear in the Orion Floe editor.
#     version = "0.1.0"
#
#     title = "SSTMAP GIST Analysis"
#
#     description = """
#         SSTMap performs Water Thermodynamics analysis.
#         SSTMaps supports hydration site analysis (HSA)
#         and Grid Inhomogeneous Solvation Theory (GIST).
#
#         SSTMap has been developed at Kurtzman Lab Lehman College
#         For more details, please visit
#         sstmap.org @ https://github.com/KurtzmanLab/SSTMap
#         """
#     classifications = [["SSTMap Analysis", "GIST"]]
#
#     tags = [tag for lists in classifications for tag in lists]
#
#     # Override defaults for some parameters
#     parameter_overrides = {
#         "memory_mb": {"default": 6000},
#         "spot_policy": {"default": "Allowed"},
#         "prefetch_count": {"default": 1},  # 1 molecule at a time
#         "item_count": {"default": 1}  # 1 molecule at a time
#     }
#
#     # The first variable passed to a parameter must always be the variable the parameter is assigned to as a string.
#     grid_res = parameter.DecimalParameter(
#         'grid_res',
#         default=0.5,
#         max_value=0.75,
#         min_value=0.2,
#         level=ADVANCED,
#         help_text='Grid resolution in A. Default: 0.5.'
#     )
#
#     grid_dim = parameter.IntegerParameter(
#         'grid_dim',
#         default=48,
#         level=ADVANCED,
#         help_text="Number of voxels in each direction. Usually grids are square, All dimensions are the same."
#                   "Default: 48."
#     )
#
#     wat_model = parameter.StringParameter(
#         'wat_model',
#         default='TIP3P',
#         choices=['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC'],
#         level=ADVANCED,
#         help_text="Water model used during the simulation. Used to set bulk density number. Default: TIP3P."
#     )
#
#     lig_res_name = parameter.StringParameter(
#         'lig_res_name',
#         default='LIG',
#         max_length=4,
#         help_text="Resname to use to identify the ligand"
#     )
#
#     start_frame = parameter.IntegerParameter(
#         'start_frame',
#         default=0,
#         min_value=0,
#         level=ADVANCED,
#         help_text="Frame index to start the SSTMap analysis. Default: 0."
#
#     )
#
#     total_frames = parameter.IntegerParameter(
#         'total_frames',
#         max_value=100000,
#         default=100,
#         level=ADVANCED,
#         help_text="Total number of frames to process during the analysis. Default: 100."
#     )
#
#     ligand_port = RecordInputPort("ligand_port", initializer=True)
#
#     def begin(self):
#         self.opt = vars(self.args)
#         self.opt['Logger'] = self.log
#
#         # Generate dictionary of water models and bulk density
#         wat_model = ['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC']
#         wat_model_bulk_density = [0.0329, 0.0332, 0.0332, 0.0329, 0.0334, 0.0333, 0.0329, 0.0333]
#         self.wat_model_density_dic = dict(zip(wat_model, wat_model_bulk_density))
#
#         for ligand in self.ligand_port:
#             self.ligand = ligand.get_value(Fields.primary_molecule)
#
#     # Records are passed to this function for processing.
#     def process(self, record, port):
#         try:
#             if port == "intake":
#
#                 # Generation options dictionary
#                 opt = dict(self.opt)
#
#                 # Assigning the bulk density bas on water model
#                 opt['rho_bulk'] = self.wat_model_density_dic[opt['wat_model']]
#
#                 if not record.has_value(Fields.primary_molecule):
#                     self.log.error("Missing molecule Primary Molecule' field")
#                     self.failure.emit(record)
#                     return
#
#                 system = record.get_value(Fields.primary_molecule)
#
#                 if not record.has_value(Fields.title):
#                     self.log.warn("Missing record Title field")
#                     system_title = system.GetTitle()[0:12]
#                 else:
#                     system_title = record.get_value(Fields.title)
#
#                 if not record.has_value(Fields.id):
#                     raise ValueError("Missing ID Field")
#
#                 sys_id = record.get_value(Fields.id)
#
#                 sys_info = system_title + '_' + str(sys_id)
#
#                 # Get the MDStageRecord list from the record
#                 if record.has_value(Fields.md_stages):
#                     mdstages = record.get_value(Fields.md_stages)
#                 else:
#                     raise ValueError("Field md_stages is missing!")
#
#                 # Get the MDStageRecord for the production stage from the MDStageRecord list
#                 # That correspond to the last member of the MDStageRecord list
#                 mdstage_prod = mdstages[-1]
#
#                 # Get the MDSystemRecord for the production stage
#                 if mdstage_prod.has_value(Fields.md_system):
#                     mdsystem_prod = mdstage_prod.get_value(Fields.md_system)
#                 else:
#                     raise ValueError("Field md_system is missing!")
#
#                 # Get the PARMED object from the MDSystemRecord
#                 prod_topology_parmed = mdsystem_prod.get_value(Fields.structure)
#
#                 with TemporaryDirectory() as output_directory:
#                     opt['Logger'].info("{} - Temporary directory: {}\n".format(sys_info , output_directory))
#
#                     # Get the name of the trajectory from the  production MDStageRecord
#                     if mdstage_prod.has_value(Fields.trajectory):
#
#                         if in_orion():
#                             prod_traj_path = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                             prod_traj_filename = os.path.join(output_directory, "trajectory.h5")
#                             copyfile(prod_traj_path, prod_traj_filename)
#                         else:
#                             prod_traj_filename = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                     else:
#                         raise ValueError("MD_stages do not have a trajectory!")
#
#                     opt['Logger'].warn("{} - Trajectory file name: {}".format(sys_info, prod_traj_filename))
#
#                     # Get the final structure from the production stage
#                     prod_coord_eomol = mdsystem_prod.get_value(Fields.topology)
#
#                     # Extract the ligand from the final frame
#                     prot, lig, wat, excp = split(prod_coord_eomol)
#
#                     self.log.info("System name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
#                                   "Water atom numbers = {}\nExcipients atom numbers = {}".format(sys_info,
#                                                                                                  prot.NumAtoms(),
#                                                                                                  lig.NumAtoms(),
#                                                                                                  wat.NumAtoms(),
#                                                                                                  excp.NumAtoms()))
#                     # Generate Variables needed for running SSTMap
#                     ligand_filename = os.path.join(output_directory, "ligand.pdb")
#
#                     opt['Logger'].info("{} - Generate files for GIST....".format(sys_info))
#
#                     # Generate files for GIST
#                     ofs = oechem.oemolostream(ligand_filename)
#
#                     pdb_flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
#
#                     ofs.SetFlavor(oechem.OEFormat_PDB, pdb_flavor)
#
#                     if lig.GetMaxAtomIdx() > 0:
#
#                         oechem.OEWriteConstMolecule(ofs, lig)
#
#                         ligand_align = False
#                     else:
#
#                         oechem.OEWriteConstMolecule(ofs, self.ligand)
#
#                         ligand_align = True
#                     ofs.close()
#
#                     if ligand_align:
#
#                         mdstage_setup = mdstages[0]
#
#                         if mdstage_setup.get_value(Fields.stage_type) == MDStageTypes.SETUP:
#                             mdsytem_setup = mdstage_setup.get_value(Fields.md_system)
#                             setup_topology = mdsytem_setup.get_value(Fields.topology)
#
#                     opt['Logger'].info("{} - Processing Trajectory".format(sys_info))
#
#                     # In order to get meaningful energy values we need to strip the ions from the trajectory
#                     # Also, we need to fit the protein to a reference frame to remove translations and
#                     # rotations.
#                     top_filename, parm_filename, prod_traj_filename, aligned_prot_oemol = sstmap_utils.\
#                         process_trajectory(prod_traj_filename, prod_topology_parmed,
#                                            opt['lig_res_name'], output_directory,
#                                            reference_topology=setup_topology)
#                     # Get number of frames
#                     total_number_frames = 0
#                     if prod_traj_filename.endswith(".h5"):
#                         for chunk in md.iterload(prod_traj_filename):
#                             total_number_frames += chunk.n_frames
#                     else:
#                         for chunk in md.iterload(prod_traj_filename, top=top_filename):
#                             total_number_frames += chunk.n_frames
#
#                     # Start SSTMap GIST analysis
#                     # Change to tmp directory to avoid data overwrite
#                     cwd = os.getcwd()
#                     os.chdir(output_directory)
#
#                     opt['Logger'].info("{} - Starting GIST Calculation".format(sys_info))
#
#                     # Initialize GIST calculation
#                     gist = sm.GridWaterAnalysis(topology_file=top_filename,
#                                                 trajectory=prod_traj_filename,
#                                                 start_frame=0,
#                                                 num_frames=total_number_frames,
#                                                 supporting_file=parm_filename,
#                                                 ligand_file=ligand_filename,
#                                                 grid_dimensions=[opt['grid_dim'], opt['grid_dim'], opt['grid_dim']],
#                                                 grid_resolution=[opt['grid_res'], opt['grid_res'], opt['grid_res']])
#
#                     # Create new record with results
#                     new_record = OERecord()
#                     new_record.set_value(Fields.primary_molecule, aligned_prot_oemol)
#                     new_record.set_value(Fields.title, system_title),
#                     new_record.set_value(Fields.id, sys_id)
#
#                     # Print System summary from GISt
#                     gist.print_system_summary()
#
#                     # Make GIST calculations
#                     gist.calculate_grid_quantities(hbonds=True)
#
#                     # Write GIST Data
#                     gist.write_data()
#
#                     # Generate constat to remove the density weight
#                     g_const = opt['grid_res'] * opt['grid_res'] * opt['grid_res']
#
#                     # Extract voxel coordinates
#                     x_gist_voxels = gist.voxeldata[:, GISTFields.x]
#                     y_gist_voxels = gist.voxeldata[:, GISTFields.y]
#                     z_gist_voxels = gist.voxeldata[:, GISTFields.z]
#
#                     # Extract the data use to generate the OEGrids
#                     g_gO = gist.voxeldata[:, GISTFields.g_O]
#                     g_gH = gist.voxeldata[:, GISTFields.g_H]
#                     g_Eww = gist.voxeldata[:, GISTFields.E_ww_dens]
#                     g_Esw = gist.voxeldata[:, GISTFields.E_sw_dens]
#                     g_So = gist.voxeldata[:, GISTFields.TS_or_dens]
#                     g_St = gist.voxeldata[:, GISTFields.TS_tr_dens]
#
#                     # Remove the density weight
#                     for g_data in [g_St, g_So, g_Eww, g_Esw]:
#                         g_data = g_data * g_const
#
#                     # Calculate total Energy = Ewat-wat + Esolute-wat
#                     g_Etot = g_Esw + g_Eww
#
#                     # Calculate total Entropy = Sorient + Strans
#                     g_Stot = g_So + g_St
#
#                     # Calculate Helmholtz free energy
#                     g_A = g_Etot - g_Stot
#
#                     # create dict for calc values
#                     calc_values = {97: 'Etot',
#                                    98: 'Stot',
#                                    99: 'FreeE'}
#
#                     # Write data to OEGrid
#                     grid_data_field_num = [GISTFields.g_O, GISTFields.g_H,
#                                            GISTFields.E_ww_dens, GISTFields.E_sw_dens,
#                                            GISTFields.TS_or_dens, GISTFields.TS_tr_dens, 97, 98, 99]
#
#                     grid_data = [g_gO, g_gH, g_Eww, g_Esw, g_So, g_St, g_Etot, g_Stot, g_A]
#
#                     grid_data_comb = list(zip(grid_data, grid_data_field_num))
#
#                     # Get grid center
#                     grid_center = gist.center.tolist()
#                     grid_orig = gist.origin.tolist()
#
#                     # Write grid information necessary to recreate the grids from the summary
#                     grid_info_fn = os.path.join(output_directory, "gist_grid_data.txt")
#
#                     with open(grid_info_fn, "w") as g_ofs:
#                         g_ofs.write("grid dimensions: {} {} {}\ngrid center: {} {} {}\ngrid origin: {} {} {}\n"
#                                     "grid resolution: {}\n".format(opt['grid_dim'], opt['grid_dim'], opt['grid_dim'],
#                                                                    grid_center[0], grid_center[1], grid_center[2],
#                                                                    grid_orig[0], grid_orig[1], grid_orig[2],
#                                                                    opt['grid_res']))
#                     for data, data_num in grid_data_comb:
#
#                         if data_num < 90:
#                             grid_name = GISTFields.data_titles[data_num]
#                         else:
#                             grid_name = calc_values[data_num]
#
#                         # Initializing the OEGrid object
#                         grid = oegrid.OEScalarGrid(opt['grid_dim'], opt['grid_dim'], opt['grid_dim'],
#                                                    grid_center[0],
#                                                    grid_center[1], grid_center[2], opt['grid_res'])
#
#                         grid.SetTitle(grid_name)
#
#                         # Set values in OEGrids
#                         for data_pnt in range(grid.GetSize()):
#                             x = x_gist_voxels[data_pnt]
#                             y = y_gist_voxels[data_pnt]
#                             z = z_gist_voxels[data_pnt]
#                             grid.SetValue(x, y, z, data[data_pnt])
#
#                         grid_field = OEField(grid_name, Types.Chem.Grid)
#                         new_record.set_value(grid_field, grid)
#
#                         # Write OEGrid
#                         grid_file_name = grid_name + ".grd"
#                         grid_file_name_wpath = os.path.join(output_directory, grid_file_name)
#                         oegrid.OEWriteGrid(grid_file_name_wpath, grid)
#
#                     # Write dx file of all calculated quantities
#                     gist.generate_dx_files()
#
#                     # Print Calculation Summary
#                     gist.print_calcs_summary()
#
#                     gist_data = os.path.join(cwd, "GIST_Results_data")
#                     shutil.copytree(output_directory, gist_data)
#
#                     self.success.emit(new_record)
#
#         except:
#             # Attach an error message to the molecule that failed
#             self.log.error(traceback.format_exc())
#             print(traceback.format_exc(), flush=True)
#             # Return failed mol
#             self.failure.emit(record)
