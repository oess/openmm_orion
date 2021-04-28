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

from MDOrion.Standards import Fields, MDStageNames

from pathlib import Path

from floereport import FloeReport, LocalFloeReport

from orionclient.session import in_orion, OrionSession, get_session

from orionclient.types import File

import os
from os import environ

import numpy as np

import MDOrion.TrjAnalysis.utils as utl

import oetrajanalysis.OETrajBasicAnalysis_utils as oetrjutl
import oetrajanalysis.Clustering_utils as clusutl

import ensemble2img

from tempfile import TemporaryDirectory

from openeye import oechem

from openeye import oedepict

import traceback

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from MDOrion.Standards.mdrecord import MDDataRecord

import MDOrion.TrjAnalysis.STMD_FloeReport_utils as flrpt

import tarfile

import shutil

from MDOrion.TrjAnalysis.dowload_fr_utils import download_report


class MDFloeReportCube(RecordPortsMixin, ComputeCube):
    title = "MDFloeReportCube"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Report']
    description = """
    The floe report cube generates an Orion floe report tiling the input ligands.
    Each input record must have ligand ID, ligand title, ligand name, the ligand
    depiction as svg string, the html report string linked to the ligand and
    optionally the ligand report label.
    """

    uuid = "58a012d2-69e9-4d15-ba17-66f65c55dec5"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    floe_report_title = parameters.StringParameter("floe_report_title", default="MD Report")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.floe_report_dic = dict()

        if in_orion():
            job_id = environ.get('ORION_JOB_ID')
            self.floe_report = FloeReport.start_report(self.opt['floe_report_title'], job_id=job_id)
        else:
            self.floe_report = LocalFloeReport.start_report(self.opt['floe_report_title'])

    def process(self, record, port):

        try:

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            system_title = mdrecord.get_title

            if mdrecord.has_conf_id:
                sort_key = (1000 * mdrecord.get_lig_id) + mdrecord.get_conf_id
            else:
                sort_key = mdrecord.get_lig_id

            if not record.has_value(Fields.floe_report):
                raise ValueError("Missing the report field for the system {}".format(system_title))

            report_string = record.get_value(Fields.floe_report)

            if not record.has_value(Fields.ligand_name):
                raise ValueError("Missing the ligand name field")

            ligand_name = record.get_value(Fields.ligand_name)
            if len(ligand_name) < 15:
                page_title = ligand_name
            else:
                page_title = ligand_name[0:13] + '...'

            if not record.has_value(Fields.floe_report_svg_lig_depiction):
                raise ValueError("Missing the ligand  depiction field")

            ligand_svg = record.get_value(Fields.floe_report_svg_lig_depiction)

            if not record.has_value(Fields.floe_report_label):
                floe_report_label = ligand_name
            else:
                floe_report_label = record.get_value(Fields.floe_report_label)

            page = self.floe_report.create_page(page_title, is_index=False)
            page_link = page.get_link()
            page.set_from_string(report_string)

            record.set_value(Fields.floe_report_URL, page_link)

            self.floe_report_dic[sort_key] = (page_link, ligand_svg, floe_report_label)

            if in_orion():
                record.set_value(Fields.floe_report_collection_id, self.floe_report.collection.id)

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

                page_link, ligand_svg, label = self.floe_report_dic[key]

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
            self.opt['Logger'].warn("It was not possible to generate the floe report: {}".format(str(e)))

        return


class ClusterOETrajCube(RecordPortsMixin, ComputeCube):
    title = 'Cluster Ligand Traj OEMol'
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Clustering', 'Ligand', 'Protein']

    description = """
    Cluster Ligand multiconf MD trajectory OEMol

    This Cube will take in the MD traj OEMols containing
    the protein and ligand components of the complex and cluster
    them based on ligand RMSD.
    """

    uuid = "b503c2f4-12e6-49c7-beb6-ee17da177ec2"

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
            opt['Logger'].info(' Beginning ClusterOETrajCube')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to cluster MD Traj'
                .format(system_title) )

            # Get the ligand which will be a multiconformer molecule with the starting
            # conformer for each simulation
            if not record.has_field(Fields.ligand):
                raise ValueError('{} could not find the ligand field'.format(system_title))
            ligand = utl.RequestOEFieldType(record, Fields.ligand)
            lig_name = utl.RequestOEFieldType(record, Fields.ligand_name)

            # Get the poseId vector which addresses each frame of the trajectory to its
            # parent starting conformer.
            if not record.has_field(Fields.Analysis.poseIdVec):
                raise ValueError('{} could not find the poseId vector'.format(system_title))
            poseIdVec = utl.RequestOEFieldType(record, Fields.Analysis.poseIdVec)

            # Get the ligand trajectory OEMol with one conformer per trajectory frame
            if not record.has_field(Fields.Analysis.oetraj_rec):
                raise ValueError('{} could not find the traj record'.format(system_title))
            opt['Logger'].info('{} found the traj record'.format(system_title))
            oetrajRecord = record.get_value(Fields.Analysis.oetraj_rec)
            ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} : got ligTraj with {} atoms, {} confs'.format(
                system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )

            # Cluster ligand trajs into a clustering results dictionary by RMSD and rotBond features
            opt['Logger'].info('{} starting clustering {} traj frames by RMSD and rotBond features'.format(
                system_title, len(poseIdVec)) )
            torScale = 0.5
            epsScal = 0.05
            clusResults = clusutl.ClusterLigTrajDBSCAN(ligand, poseIdVec, ligTraj, torScale, epsScal)

            opt['Logger'].info('{} clustering completed finding {} clusters with {} outliers'.format(
                system_title, clusResults['nClusters'], clusResults['nOutliers']) )
            opt['Logger'].info('{} cluster counts: {}'.format(
                system_title, clusResults['ClusterCounts']) )

            # Clusters are ordered by decreasing size, so first n clusters will be major clusters
            clusterCounts = clusResults['ClusterCounts']
            majorClusThreshold = 0.05
            clusResults['MajorClusThreshold'] = majorClusThreshold
            nLarge = 0
            for count in clusterCounts:
                if count/clusResults['nFrames'] >= majorClusThreshold:
                    nLarge += 1
            clusResults['nMajorClusters'] = nLarge

            # Create new record with trajClus results
            trajClus = OERecord()
            #
            # store trajClus results dict (Plain Old Data only) on the record as a JSON object
            trajClus.set_value(Fields.Analysis.oeclus_dict, clusResults)
            opt['Logger'].info('{} Saved clustering results in dict with keys:'.format(system_title) )
            for key in clusResults.keys():
                opt['Logger'].info('{} : TrajClusDict key {}'.format(system_title, key) )

            # Set the TrajClus record on the top-level record
            record.set_value(Fields.Analysis.oeclus_rec, trajClus)
            opt['Logger'].info('{} cluster results written to TrajClus OERecord'.format(system_title) )

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ClusterOETrajCube on {}'.format(str(e), system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class MakeClusterTrajOEMols(RecordPortsMixin, ComputeCube):
    title = 'Make Cluster Protein and Ligand average and median OEMols'
    # version = "0.2.0"
    classification = [["Analysis"]]
    tags = ['Clustering', 'Ligand', 'Protein']

    description = """
    Make multiconf MD trajectory OEMols for Ligand and Protein by Cluster

    This Cube will use the clustering results in conjunction with the
    MD trajectory OEMols for protein and ligand to generate per-cluster
    protein and ligand average and median structures.
    """

    uuid = "d500eaf0-b775-4caf-af5e-a6d79f5c0603"

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
            opt['Logger'].info(' Beginning MakeClusterTrajOEMols Cube')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to make trajectory OEMols for Ligand and Protein by Cluster'
                               .format(system_title))

            if not record.has_field(Fields.Analysis.oeclus_rec):
                raise ValueError('{} does not have TrajClus record'.format(system_title))
            else:
                opt['Logger'].info('{} found TrajClus record'.format(system_title))

            trajClusRecord = utl.RequestOEFieldType(record, Fields.Analysis.oeclus_rec)

            # Get the cluster info dict off the oeclus_dict field
            if not trajClusRecord.has_field(Fields.Analysis.oeclus_dict):
                raise ValueError('{} could not find the oeclus_dict field'.format(system_title))
            else:
                opt['Logger'].info('{} found the oeclus_dict field'.format(system_title))

            # Extract the relevant clustering information from the trajClus results dict
            trajClus = trajClusRecord.get_value(Fields.Analysis.oeclus_dict)
            opt['Logger'].info('{} retrieved Cluster info on {} frames giving {} clusters'
                               .format(system_title, trajClus['nFrames'], len(trajClus['ClusterCounts'])))

            # Extract the ligTraj OEMol from the OETraj record
            if not record.has_field(Fields.Analysis.oetraj_rec):
                raise ValueError('{} could not find the traj record'.format(system_title))
            else:
                opt['Logger'].info('{} found the traj record'.format(system_title))
            oetrajRecord = record.get_value(Fields.Analysis.oetraj_rec)
            ligTraj = utl.RequestOEField(oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} got ligTraj with {} atoms, {} confs'.format(
                system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()))

            # Extract the protTraj OEMol from the OETraj record
            mdtrajrecord = MDDataRecord(oetrajRecord)
            protTraj = mdtrajrecord.get_protein_traj
            opt['Logger'].info('{} got protTraj with {} atoms, {} confs'.format(
                system_title, protTraj.NumAtoms(), protTraj.NumConfs()) )
            del mdtrajrecord

            # Generate average and median protein and ligand OEMols from ligTraj, protTraj
            opt['Logger'].info('{} Generating entire trajectory median and average OEMols for protein and ligand '.
                               format(system_title))
            ligMedian, protMedian, ligAverage, protAverage = oetrjutl.AnalyseProteinLigandTrajectoryOEMols(
                ligTraj, protTraj)

            # Add prot and lig medians and averages to OETraj record
            #oetrajRecord.set_value(OEField('LigMedian', Types.Chem.Mol), ligMedian)
            #oetrajRecord.set_value(OEField('ProtMedian', Types.Chem.Mol), protMedian)
            oetrajRecord.set_value(OEField('LigAverage', Types.Chem.Mol), ligAverage)
            #oetrajRecord.set_value(OEField('ProtAverage', Types.Chem.Mol), protAverage)

            # Generate interactive trajectory SVG for the whole trajectory and place on oetrajRecord
            opt['Logger'].info('{} Generating entire trajectory interactive SVG'.format(system_title))
            trajSVG = ensemble2img.run_ensemble2img(ligMedian, protMedian, ligTraj, protTraj)
            TrajSVG_field = OEField('TrajSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            oetrajRecord.set_value(TrajSVG_field, trajSVG)
            record.set_value(Fields.Analysis.oetraj_rec, oetrajRecord)

            # Generate per-cluster info for major clusters
            nMajorClusters = trajClus['nMajorClusters']
            byClusTrajSVG = []
            if nMajorClusters > 0:
                opt['Logger'].info('{}: Making cluster mols for {} major clusters'.format(system_title,nMajorClusters))
                clusLigAvgMol = oechem.OEMol(ligTraj)
                clusLigAvgMol.DeleteConfs()
                clusProtAvgMol = oechem.OEMol(protTraj)
                clusProtAvgMol.DeleteConfs()
                clusLigMedMol = oechem.OEMol(ligTraj)
                clusLigMedMol.DeleteConfs()
                clusProtMedMol = oechem.OEMol(protTraj)
                clusProtMedMol.DeleteConfs()

                # for each major cluster generate SVG and average and median for protein and ligand
                for clusID in range(nMajorClusters):
                    opt['Logger'].info('Extracting cluster {} from {}'.format( clusID, system_title ))
                    clusLig = clusutl.TrajOEMolFromCluster( ligTraj, trajClus['ClusterVec'], clusID)
                    opt['Logger'].info( 'ligand cluster {} with {} confs'.format(clusID,clusLig.NumConfs()) )
                    clusProt = clusutl.TrajOEMolFromCluster( protTraj, trajClus['ClusterVec'], clusID)
                    opt['Logger'].info('protein cluster {} with {} confs'.format(clusID,clusProt.NumConfs()) )
                    opt['Logger'].info('generating representative protein average and median confs')
                    #
                    ligMed, protMed, ligAvg, protAvg = oetrjutl.AnalyseProteinLigandTrajectoryOEMols( clusLig, clusProt)
                    confTitle = 'clus '+str(clusID)
                    conf = clusLigAvgMol.NewConf(ligAvg)
                    conf.SetTitle(confTitle)
                    conf = clusProtAvgMol.NewConf(protAvg)
                    conf.SetTitle(confTitle)
                    conf = clusLigMedMol.NewConf(ligMed)
                    conf.SetTitle(confTitle)
                    conf = clusProtMedMol.NewConf(protMed)
                    conf.SetTitle(confTitle)
                    #
                    opt['Logger'].info('generating cluster SVG for cluster {}'.format(clusID) )
                    clusSVG = ensemble2img.run_ensemble2img(ligAvg, protAvg, clusLig, clusProt)
                    byClusTrajSVG.append(clusSVG)

                # style the molecules and put the results on trajClus record
                clusLigAvgMol.SetTitle('Average '+clusLigAvgMol.GetTitle())
                clusProtAvgMol.SetTitle('Average '+clusProtAvgMol.GetTitle())
                utl.StyleTrajProteinLigandClusters(clusProtAvgMol,clusLigAvgMol)
                trajClusRecord.set_value(Fields.Analysis.ClusLigAvg_fld, clusLigAvgMol)
                trajClusRecord.set_value(Fields.Analysis.ClusProtAvg_fld, clusProtAvgMol)
                #
                clusLigMedMol.SetTitle('Median '+clusLigMedMol.GetTitle())
                clusProtMedMol.SetTitle('Median '+clusProtMedMol.GetTitle())
                utl.StyleTrajProteinLigandClusters(clusProtMedMol,clusLigMedMol)
                trajClusRecord.set_value(Fields.Analysis.ClusLigMed_fld, clusLigMedMol)
                trajClusRecord.set_value(Fields.Analysis.ClusProtMed_fld, clusProtMedMol)

            # case when no major clusters are found
            else:  # number of clusters is zero
                opt['Logger'].info('No major clusters found for {}'.format(system_title))
                # In lieu of cluster mols, copy traj and median mols to top level
                #clusLigAvgMol = utl.RequestOEField( oetrajRecord, 'LigAverage', Types.Chem.Mol)
                #clusProtAvgMol = utl.RequestOEField( oetrajRecord, 'ProtAverage', Types.Chem.Mol)
                clusLigAvgMol = ligAverage
                clusProtAvgMol = protAverage
                clusLigAvgMol.SetTitle('Average '+clusLigAvgMol.GetTitle())
                clusProtAvgMol.SetTitle('Average '+clusProtAvgMol.GetTitle())
                utl.SetProteinLigandVizStyle(clusProtAvgMol, clusLigAvgMol)
                #
                #clusLigMedMol = utl.RequestOEField( oetrajRecord, 'LigMedian', Types.Chem.Mol)
                #clusProtMedMol = utl.RequestOEField( oetrajRecord, 'ProtMedian', Types.Chem.Mol)
                clusLigMedMol = ligMedian
                clusProtMedMol = protMedian
                clusLigMedMol.SetTitle('Median '+clusLigMedMol.GetTitle())
                clusProtMedMol.SetTitle('Median '+clusProtMedMol.GetTitle())
                utl.SetProteinLigandVizStyle(clusProtMedMol, clusLigMedMol)

            # Set prot and lig clus average mols on top-level record for 3D vis
            record.set_value(Fields.Analysis.ClusLigAvg_fld, clusLigAvgMol)
            record.set_value(Fields.Analysis.ClusProtAvg_fld, clusProtAvgMol)
            record.set_value(Fields.Analysis.ClusLigMed_fld, clusLigMedMol)
            record.set_value(Fields.Analysis.ClusProtMed_fld, clusProtMedMol)

            # put highlighted carbon styling on the default molecule (initial pose ligand)
            primaryMol = utl.RequestOEFieldType( record, Fields.primary_molecule)
            utl.HighlightStyleMolecule(primaryMol)
            record.set_value(Fields.primary_molecule, primaryMol)

            byClusTrajSVG_field = OEField('ByClusTrajSVG', Types.StringVec)
            trajClusRecord.set_value(byClusTrajSVG_field, byClusTrajSVG)

            # Set the TrajClus record on the top-level record
            record.set_value(Fields.Analysis.oeclus_rec, trajClusRecord)
            opt['Logger'].info('{} finished adding to trajClus OERecord'.format(system_title) )

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ClusterOETrajCube on {}'.format(str(e), system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class ClusterPopAnalysis(RecordPortsMixin, ComputeCube):
    title = 'Population Analysis of Traj Clusters'
    # version = "0.1.0"
    classification = [["Analysis"]]
    tags = ['Clustering', 'Ligand', 'Trajectory']

    description = """
    Population Analysis of Traj Clusters

    This Cube analyzes the ligand trajectory cluster in terms of
    their occurrence, size, and proximity to the starting pose(s).
    """

    uuid = "ddffa0d3-6958-4df6-b322-5020e049a040"

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
            opt['Logger'].info(' Beginning ClusterPopAnalysis')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting Population Analysis of Traj Clusters'
                .format(system_title) )

            # Get the ligand which will be a multiconformer molecule comprising the parent
            # starting conformer for each separate MD simulation
            if not record.has_field(Fields.ligand):
                raise ValueError('{} could not find the ligand field'.format(system_title))
            ligand = utl.RequestOEFieldType(record, Fields.ligand)
            lig_name = utl.RequestOEFieldType(record, Fields.ligand_name)

            # Extract the ligTraj OEMol from the OETraj record
            if not record.has_field(Fields.Analysis.oetraj_rec):
                raise ValueError('{} could not find the traj record'.format(system_title))
            oetrajRecord = record.get_value(Fields.Analysis.oetraj_rec)
            ligTraj = utl.RequestOEField(oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} got ligTraj with {} atoms, {} confs'.format(
                system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()))

            # Get the poseId vector which addresses each frame of the trajectory to its
            # parent starting conformer.
            if not record.has_field(Fields.Analysis.poseIdVec):
                raise ValueError('{} could not find the poseId vector'.format(system_title))
            poseIdVec = utl.RequestOEFieldType(record, Fields.Analysis.poseIdVec)

            # Get the clustering results dict from the traj clustering record
            if not record.has_field(Fields.Analysis.oeclus_rec):
                raise ValueError('{} could not find the cluster record'.format(system_title))
            opt['Logger'].info('{} found the cluster record'.format(system_title))
            oeclusRecord = record.get_value(Fields.Analysis.oeclus_rec)
            clusResults = utl.RequestOEFieldType( oeclusRecord, Fields.Analysis.oeclus_dict)
            #for key in clusResults.keys():
            #    opt['Logger'].info('{} : clusResults key {}'.format(system_title, key) )

            # Generate the fractional cluster populations by conformer, and conformer populations by cluster
            popResults = clusutl.AnalyzeClustersByConfs(ligand, poseIdVec, clusResults)

            # Get the PBSA data dict from the record
            if not record.has_field(Fields.Analysis.oepbsa_dict):
                raise ValueError('{} could not find the PBSA JSON object'.format(system_title))
            opt['Logger'].info('{} found the PBSA JSON record'.format(system_title))
            PBSAdata = utl.RequestOEFieldType(record,Fields.Analysis.oepbsa_dict)
            #for key in PBSAdata.keys():
            #    opt['Logger'].info('{} : PBSAdata key {} {}'.format(system_title, key, len(PBSAdata[key])) )

            # Generate the cluster MMPBSA mean and standard error
            MMPBSAbyClus = clusutl.MeanSerrByClusterEnsemble(popResults, PBSAdata['OEZap_MMPBSA6_Bind'])
            popResults['OEZap_MMPBSA6_ByClusMean'] = MMPBSAbyClus['ByClusMean']
            popResults['OEZap_MMPBSA6_ByClusSerr'] = MMPBSAbyClus['ByClusSerr']
            popResults['OEZap_MMPBSA6_ByConfMean'] = MMPBSAbyClus['ByConfMean']
            popResults['OEZap_MMPBSA6_ByConfSerr'] = MMPBSAbyClus['ByConfSerr']

            # Get the TrajBintScore data vector from the BintScore record
            if not record.has_field(Fields.Bint.oebint_rec):
                raise ValueError('{} could not find the BintScore record'.format(system_title))
            opt['Logger'].info('{} found the BintScore record'.format(system_title))
            oebint_rec = record.get_value(Fields.Bint.oebint_rec)
            trajBintScoreList = utl.RequestOEFieldType(oebint_rec, Fields.Bint.trajBintScoreList)

            # Generate the cluster TrajBintScore mean and standard error
            trajBintScorebyClus = clusutl.MeanSerrByClusterEnsemble(popResults, trajBintScoreList)
            popResults['TrajBintScore_ByClusMean'] = trajBintScorebyClus['ByClusMean']
            popResults['TrajBintScore_ByClusSerr'] = trajBintScorebyClus['ByClusSerr']
            popResults['TrajBintScore_ByConfMean'] = trajBintScorebyClus['ByConfMean']
            popResults['TrajBintScore_ByConfSerr'] = trajBintScorebyClus['ByConfSerr']

            # Generate by-cluster mean and serr RMSDs to the starting confs
            ClusRMSDByConf = clusutl.ClusterRMSDByConf(ligand, ligTraj, clusResults)
            popResults['confRMSDsByClusMean'] = ClusRMSDByConf['confRMSDsByClusMean']
            popResults['confRMSDsByClusSerr'] = ClusRMSDByConf['confRMSDsByClusSerr']

            # Put these results on the record as a POD JSON object
            oeclusRecord.set_value(Fields.Analysis.cluspop_dict, popResults)
            record.set_value(Fields.Analysis.oeclus_rec, oeclusRecord)

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ClusterPopAnalysis on {}'.format(str(e), system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class TrajAnalysisReportDataset(RecordPortsMixin, ComputeCube):
    title = 'Prepare the Traj Analysis Dataset for Reporting'
    # version = "0.1.0"
    classification = [["Analysis"]]
    tags = ['Clustering', 'Ligand', 'Trajectory']

    description = """
    Prepare the Traj Analysis Dataset for Reporting

    This reorganizes and modifies the Short Traj MD Analysis datarecord to
    trim down and focus the results for display in Orion the associated Floe Report.
    """

    uuid = "f4e2359b-aab5-499c-a4fc-735856fea93b"

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
            opt['Logger'].info(' Beginning TrajAnalysisReportDataset')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Prepare the Traj Analysis Dataset for Report'
                               .format(system_title))

            # Get the ligand which will be a multiconformer molecule with the starting
            # conformer for each simulation
            if not record.has_field(Fields.ligand):
                raise ValueError('{} could not find the ligand field'.format(system_title))
            ligand = utl.RequestOEFieldType(record, Fields.ligand)
            lig_name = utl.RequestOEFieldType(record, Fields.ligand_name)

            # Get the poseId vector which addresses each frame of the trajectory to its
            # parent starting conformer.
            if not record.has_field(Fields.Analysis.poseIdVec):
                raise ValueError('{} could not find the poseId vector'.format(system_title))
            poseIdVec = utl.RequestOEFieldType(record, Fields.Analysis.poseIdVec)

            # Copy InitBintScore and TrajBintScore from the oebint_rec record to the top level
            if not record.has_field(Fields.Bint.oebint_rec):
                raise ValueError('{} could not find the BintScore record'.format(system_title))
            opt['Logger'].info('{} found the BintScore record'.format(system_title))
            bintRecord = record.get_value(Fields.Bint.oebint_rec)
            initBintScore = utl.RequestOEFieldType(bintRecord, Fields.Bint.initBintScore)
            trajBintScore = utl.RequestOEFieldType(bintRecord, Fields.Bint.trajBintScore)
            trajBintStderr = utl.RequestOEFieldType(bintRecord, Fields.Bint.trajBintStderr)
            record.set_value(Fields.Bint.initBintScore, initBintScore)
            record.set_value(Fields.Bint.trajBintScore, trajBintScore)
            record.set_value(Fields.Bint.trajBintStderr, trajBintStderr)

            # Get the clustering results dict from the traj clustering record
            if not record.has_field(Fields.Analysis.oeclus_rec):
                raise ValueError('{} could not find the cluster record'.format(system_title))
            opt['Logger'].info('{} found the cluster record'.format(system_title))
            oeclusRecord = record.get_value(Fields.Analysis.oeclus_rec)
            clusResults = utl.RequestOEFieldType(oeclusRecord, Fields.Analysis.oeclus_dict)
            #for key in clusResults.keys():
            #    opt['Logger'].info('{} : clusResults key {}'.format(system_title, key))

            # Set the number of major clusters at the top level
            nMajorClusters = clusResults['nMajorClusters']
            record.set_value(Fields.Analysis.n_major_clusters, nMajorClusters)

            # Generate simple plots for floe report
            opt['Logger'].info('{} plotting cluster strip plot'.format(system_title) )
            trajClus_svg = clusutl.ClusterMembersshipPlot(clusResults, poseIdVec)
            #trajClus_svg = clusutl.ClusterLigTrajClusPlot(clusResults)
            ClusSVG_field = OEField( 'TrajClusSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            oeclusRecord.set_value( ClusSVG_field, trajClus_svg)

            # Calculate RMSD of ligand traj from ligand initial pose
            #ligInitPose = utl.RequestOEFieldType(record, Fields.ligand)
            #vecRmsd = oechem.OEDoubleArray(ligTraj.GetMaxConfIdx())
            #oechem.OERMSD(ligInitPose, ligTraj, vecRmsd)
            #trajClus.set_value(Fields.Analysis.lig_traj_rmsd, list(vecRmsd) )
            #opt['Logger'].info('{} plotting strip plot of ligand RMSD from initial pose'.format(system_title) )
            #rmsdInit_svg = clusutl.RmsdFromInitialPosePlot( clusResults['ClusterVec'], vecRmsd)
            # Put simple plot results on trajClus record
            #
            #rmsdInit_field = OEField( 'rmsdInitPose', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            #trajClus.set_value(rmsdInit_field, rmsdInit_svg)
            #

            # Set the TrajClus record on the top-level record
            record.set_value(Fields.Analysis.oeclus_rec, oeclusRecord)
            opt['Logger'].info('{} added report info to oeclusRecord OERecord'.format(system_title) )

            # This last section is about setting the trajectory ensemble MMPBSA value. There are 2 cases:
            # 1) There is at least one major cluster, so make a Boltzmann-weighted average of all major clusters
            # 2) There are no major clusters, so make a simple ensemble average of the whole traj.
            floe_report_label = ""
            #
            if nMajorClusters>0:
                # 1) There is at least one major cluster, so make a Boltzmann-weighted average of all major clusters:
                #
                # Get the results dict for the Cluster Population analysis
                if not oeclusRecord.has_field(Fields.Analysis.cluspop_dict):
                    raise ValueError('{} could not find the clusConf population JSON object'.format(system_title))
                opt['Logger'].info('{} found the clusConf population JSON record'.format(system_title))
                opt['Logger'].info('{} calculating Boltzmann-weighted MMPBSA average'.format(system_title))
                popResults = utl.RequestOEFieldType(oeclusRecord, Fields.Analysis.cluspop_dict)
                if 'OEZap_MMPBSA6_ByClusMean' not in popResults.keys():
                    raise ValueError('{} could not find OEZap_MMPBSA6_ByClusMean in popResults'.format(system_title))
                #
                # If >1 field, exclude the last field in the energy vector since it is Outliers+MinorClusters
                if len(popResults['OEZap_MMPBSA6_ByClusMean'])>1:
                    MMPBSA6_ByClusMean = np.array(popResults['OEZap_MMPBSA6_ByClusMean'][:-1])
                    MMPBSA6_ByClusSerr = np.array(popResults['OEZap_MMPBSA6_ByClusSerr'][:-1])
                else:
                    MMPBSA6_ByClusMean = np.array(popResults['OEZap_MMPBSA6_ByClusMean'])
                    MMPBSA6_ByClusSerr = np.array(popResults['OEZap_MMPBSA6_ByClusSerr'])
                # Calculate the Boltzmann-weighted probabilities for each cluster as a numpy array
                BoltzWtProb = oetrjutl.BoltzmannWeightedProbabilities(MMPBSA6_ByClusMean)
                boltzMean = (MMPBSA6_ByClusMean * BoltzWtProb).sum()
                boltzSerr = (MMPBSA6_ByClusSerr * BoltzWtProb).sum()
                # Add to the record the MMPBSA mean and std
                record.set_value(Fields.Analysis.mmpbsa_traj_mean, boltzMean)
                record.set_value(Fields.Analysis.mmpbsa_traj_serr, boltzSerr)
                # Add to the record the Average MMPBSA floe report label
                floe_report_label = "MMPBSA score:<br>{:.1f}  &plusmn; {:.1f} kcal/mol".format(
                    boltzMean, boltzSerr)

            else:
                # 2) There are no major clusters, so make a simple ensemble average of the whole traj.
                # Get the PBSA data dict from the record
                if not record.has_field(Fields.Analysis.oepbsa_dict):
                    raise ValueError('{} could not find the PBSA data JSON object'.format(system_title))
                opt['Logger'].info('{} found the PBSA data JSON record'.format(system_title))
                opt['Logger'].info('{} calculating simple ensemble MMPBSA average'.format(system_title))
                PBSAdata = utl.RequestOEFieldType(record, Fields.Analysis.oepbsa_dict)
                if 'OEZap_MMPBSA6_Bind' not in PBSAdata.keys():
                    raise ValueError('{} could not find OEZap_MMPBSA6_Bind in PBSAdata'.format(system_title))

                # Clean MMPBSA mean and serr to avoid nans and high zap energy values
                avg_mmpbsa, serr_mmpbsa = clusutl.clean_mean_serr(PBSAdata['OEZap_MMPBSA6_Bind'])

                # Add to the record the MMPBSA mean and std
                record.set_value(Fields.Analysis.mmpbsa_traj_mean, avg_mmpbsa)
                record.set_value(Fields.Analysis.mmpbsa_traj_serr, serr_mmpbsa)
                # Add to the record the Average MMPBSA floe report label
                floe_report_label = "MMPBSA score:<br>{:.1f}  &plusmn; {:.1f} kcal/mol".format(
                    avg_mmpbsa, serr_mmpbsa)

            # Revise floe report label to include nMajorClusters
            floe_report_label = "# clusters: " + str(nMajorClusters) + "<br>" + floe_report_label
            record.set_value(Fields.floe_report_label, floe_report_label)

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ClusterPopAnalysis on {}'.format(str(e), system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class MDTrajAnalysisClusterReport(RecordPortsMixin, ComputeCube):
    title = 'Extract relevant outputs of MD Traj Cluster  Analysis'
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Ligand', 'Protein']

    description = """
    Extract relevant outputs of Ligand and Protein
    Short Traj MD Traj Analysis and write them to files.

    This Cube takes as input the OERecord containing the work
    product of trajectory analysis on Short Traj MD results.
    """

    uuid = "42f2eef0-60aa-46f8-8d55-c8f10576e319"

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

            # title of entire solvated protein-ligand system
            opt['Logger'].info('Starting Floe Report generation for MD Traj Analysis')

            system_title = utl.RequestOEFieldType(record, Fields.title)

            opt['Logger'].info('{} Attempting to extract MD Traj Analysis results'.format(system_title))

            ligInitPose = utl.RequestOEFieldType(record, Fields.ligand)

            lig_name = utl.RequestOEFieldType(record, Fields.ligand_name)

            # protInitPose = utl.RequestOEFieldType(record, Fields.protein)

            if not record.has_field(Fields.md_components):
                raise ValueError("Missing MD Components field")

            md_components = record.get_value(Fields.md_components)

            asiteSVG = utl.PoseInteractionsSVG(md_components, width=400, height=265)

            # Extract the traj SVG and Ligand average Bfactor from the OETraj record
            if not record.has_field(Fields.Analysis.oetraj_rec):
                raise ValueError('{} does not have ligand OETraj record'.format(system_title) )
            oetrajRecord = utl.RequestOEFieldType(record, Fields.Analysis.oetraj_rec)
            opt['Logger'].info('{} found OETraj record'.format(system_title))
            trajSVG = utl.RequestOEField(oetrajRecord, 'TrajSVG', Types.String)
            ligand_bfactor = utl.RequestOEField(oetrajRecord, 'LigAverage', Types.Chem.Mol)

            # Extract the three plots from the TrajClus record
            if not record.has_field(Fields.Analysis.oeclus_rec):
                raise ValueError('{} does not have TrajClus record'.format(system_title))
            # Extract the relevant traj SVG from the TrajClus record
            clusRecord = utl.RequestOEFieldType(record, Fields.Analysis.oeclus_rec)
            opt['Logger'].info('{} found TrajClus record'.format(system_title))
            trajClus_svg = utl.RequestOEField(clusRecord, 'TrajClusSVG', Types.String)
            #rmsdInit_svg = utl.RequestOEField(clusRecord, 'rmsdInitPose', Types.String)
            byClusTrajSVG = utl.RequestOEField(clusRecord, 'ByClusTrajSVG', Types.StringVec)
            opt['Logger'].info('{} found the TrajClus plots'.format(system_title))

            # Get the Clustering information
            if not clusRecord.has_field(Fields.Analysis.oeclus_dict):
                raise ValueError('{} could not find the oeclus_dict field'.format(system_title))
            else:
                opt['Logger'].info('{} found the oeclus_dict field'.format(system_title))
            # Extract the relevant clustering information from the trajClus results dict
            clusData = clusRecord.get_value(Fields.Analysis.oeclus_dict)
            opt['Logger'].info('{} found the cluster info'.format(system_title))

            # Extract the MMPBSA score for the whole trajectory
            mmpbsa_traj_mean = 999.
            mmpbsa_traj_std = 999.
            if record.has_value(Fields.Analysis.mmpbsa_traj_mean):
                mmpbsa_traj_mean = record.get_value(Fields.Analysis.mmpbsa_traj_mean)
                mmpbsa_traj_std = record.get_value(Fields.Analysis.mmpbsa_traj_serr)

            # Extract the traj BintScore for the whole trajectory
            trajbintscore = 999.
            trajbintstderr = 999.
            if record.has_value(Fields.Bint.trajBintScore):
                trajbintscore = record.get_value(Fields.Bint.trajBintScore)
                trajbintstderr = record.get_value(Fields.Bint.trajBintStderr)

            # Get the results dict for the Cluster Population analysis
            if not clusRecord.has_field(Fields.Analysis.cluspop_dict):
                raise ValueError('{} could not find the clusConf population JSON object'.format(system_title))
            opt['Logger'].info('{} found the clusConf population JSON record'.format(system_title))
            popResults = utl.RequestOEFieldType(clusRecord, Fields.Analysis.cluspop_dict)
            popTableStyles, popTableBody = flrpt.HtmlMakeClusterPopTables(popResults)

            # Make a copy of the ligand starting pose.
            # OE Prepare Depiction is removing hydrogens
            ligand_init = oechem.OEMol(ligInitPose)

            # prepare the 2D structure depiction
            oedepict.OEPrepareDepiction(ligInitPose)
            img = oedepict.OEImage(400, 300)
            oedepict.OERenderMolecule(img, ligInitPose)

            # get the palette of graph marker colors
            nClustersP1 = clusData['nClusters']+1
            clusRGB = clusutl.ColorblindRGBMarkerColors(nClustersP1)
            clusRGB[-1] = (76, 76, 76)

            with TemporaryDirectory() as output_directory:

                # write the report
                reportFName = os.path.join(output_directory, system_title + '_ClusReport.html')

                report_file = open(reportFName, 'w')

                report_file.write(flrpt._clus_floe_report_header)

                for i in range(len(byClusTrajSVG)+2):
                    report_file.write("""
                  div.cb-floe-report__tab-wrapper input:nth-of-type({clusID}):checked ~ .cb-floe-report__tab-content:nth-of-type({clusID}) {{ display: block; }}
                """.format(clusID=i+1))

                report_file.write(popTableStyles)

                report_file.write(flrpt._clus_floe_report_header2)

                report_file.write(flrpt._clus_floe_report_midHtml0.format(
                    query_depiction=oedepict.OEWriteImageToString("svg", img).decode("utf-8")))

                report_file.write("""      <h4 style="text-align: center; width: 100%">""")
                #
                # If there was no mmpbsa_traj_mean, magic number 999. so write lig_name instead
                if mmpbsa_traj_mean==999.:
                    report_file.write("""   {}</h4>""".format(lig_name))
                # Usual case: write mmpbsa_traj_mean and stderr
                else:
                    report_file.write("&lt;BintScore&gt;: {:.1f}  &plusmn; {:.1f}<br>".format(
                        trajbintscore, trajbintstderr))
                    if clusData['nClusters']>1:
                        report_file.write("&lt;MMPBSA&gt;<sup>*</sup>: {:.1f}  &plusmn; {:.1f} kcal/mol".format(
                            mmpbsa_traj_mean, mmpbsa_traj_std))
                        report_file.write("<br><sup>*</sup>Boltzmann-weighted by cluster""")
                    else:
                        report_file.write("&lt;MMPBSA&gt;: {:.1f}  &plusmn; {:.1f} kcal/mol".format(
                            mmpbsa_traj_mean, mmpbsa_traj_std))
                    report_file.write("""      </h4>""")


                analysis_txt = flrpt.MakeClusterInfoText(clusData, popResults, clusRGB)
                report_file.write("".join(analysis_txt))

                report_file.write(flrpt._clus_floe_report_midHtml1)

                report_file.write("""      <input type="radio" name="tab" id="cb-floe-report__tab-1-header" checked>
                      <label class="cb-floe-report__tab-label" for="cb-floe-report__tab-1-header">Overall</label>""")

                CurrentTabId = 1

                for i, (clus, rgb) in enumerate(zip(byClusTrajSVG, clusRGB)):
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
                      </div>""".format(traj=flrpt.trim_svg(trajSVG)))

                for clusSVG in byClusTrajSVG:
                    report_file.write("""      <div class="cb-floe-report__tab-content">
                        {traj}
                      </div>
                      """.format(traj=flrpt.trim_svg(clusSVG)))

                report_file.write("""      <div class="cb-floe-report__tab-content">
                        {traj}
                      </div>
                      """.format(traj=flrpt.trim_svg(asiteSVG)))

                #report_file.write(flrpt._clus_floe_report_midHtml2)

                report_file.write(flrpt._clus_floe_report_midHtml2a)
                report_file.write(flrpt._clus_floe_report_stripPlots.format(
                    clusters=flrpt.trim_svg(trajClus_svg)))
                    #rmsdInit=flrpt.trim_svg(rmsdInit_svg)))
                report_file.write('</div><hr style="max-width: 1000px;">')

                report_file.write(popTableBody)

                report_file.write(flrpt._clus_floe_report_Trailer)

                report_file.close()

                with open(reportFName, 'r') as f:
                    report_html_str = f.read()

                record.set_value(Fields.floe_report, report_html_str)

                # Copy Bfactors from the average Bfactor ligand to a copy of the ligand initial pose
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

                # TODO C. Bayly 2019 jul 9
                # Having written the analysis report, we know we are finished with this molecule
                # so set up the top-level record for display in Orion
                # record.set_value(Fields.primary_molecule, record.get_value(Fields.ligand))

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class ExtractMDDataCube(RecordPortsMixin, ComputeCube):
    title = "MDExtractDataCube"
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Report']
    description = """
    The cube extracts the relevant MD data generated from the Short
    Trajectory MD with Analysis floe saving the results as file in Amazon S3
    ready to be downloaded. The extracted MD data includes the protein,
    ligand and binding site water oemol trajectories, the recorded average
    and median cluster poses for the protein and ligand and the generated
    floe report.
    """

    uuid = "35a8a2b8-b765-4a6f-8d8e-6e2198feea22"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    out_file_name = parameters.StringParameter(
        'out_file_name',
        default="md_data.tar.gz",
        help_text=""""Output File name""")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.dirs = []
        self.protein_name = None
        self.floe_report_collection_id = None

    def process(self, record, port):

        try:

            # TODO CHANGE TO THE SCRATCH DIRS FOR ORION
            if not record.has_field(Fields.title):
                raise ValueError("The record is missing the title field")

            title = record.get_value(Fields.title)

            if not record.has_field(Fields.protein_name):
                raise ValueError("The record is missing the protein name field")

            if self.protein_name is None:
                self.protein_name = record.get_value(Fields.protein_name)

            if self.floe_report_collection_id is None:
                if in_orion():

                    if not record.has_field(Fields.floe_report_collection_id):
                        raise ValueError("Floe report collection id filed is missing from tje record")

                    self.floe_report_collection_id = record.get_value(Fields.floe_report_collection_id)

            Path(title).mkdir(parents=True, exist_ok=False)

            self.dirs.append(title)

            if not record.has_field(Fields.Analysis.oeclus_rec):
                raise ValueError("Cluster record field is missing")

            oeclus_rec = record.get_value(Fields.Analysis.oeclus_rec)

            for fd in oeclus_rec.get_fields():

                name = fd.get_name()

                if name in "ClusLigAvgMol":
                    clus_lig_avg = oeclus_rec.get_value(fd)
                    clus_lig_avg.SetTitle("ClusLigAvg")
                    fn = os.path.join(title, "cluster_ligand_average.oeb")
                    with oechem.oemolostream(fn) as ofs:
                        oechem.OEWriteConstMolecule(ofs, clus_lig_avg)

                elif name in "ClusLigMedMol":
                    clus_lig_med = oeclus_rec.get_value(fd)
                    clus_lig_med.SetTitle("ClusLigMed")
                    fn = os.path.join(title, "cluster_ligand_median.oeb")
                    with oechem.oemolostream(fn) as ofs:
                        oechem.OEWriteConstMolecule(ofs, clus_lig_med)

                elif name in "ClusProtAvgMol":
                    clus_prot_avg = oeclus_rec.get_value(fd)
                    clus_prot_avg.SetTitle("ClusProtAvg")
                    fn = os.path.join(title, "cluster_protein_average.oeb")
                    with oechem.oemolostream(fn) as ofs:
                        oechem.OEWriteConstMolecule(ofs, clus_prot_avg)

                elif name in "ClusProtMedMol":
                    clus_prot_med = oeclus_rec.get_value(fd)
                    clus_prot_med.SetTitle("ClusProtMed")
                    fn = os.path.join(title, "cluster_protein_median.oeb")
                    with oechem.oemolostream(fn) as ofs:
                        oechem.OEWriteConstMolecule(ofs, clus_prot_med)

            if not record.has_field(Fields.Analysis.oetraj_rec):
                raise ValueError("Multi Conf Trajectory record field is missing")

            oetraj_rec = record.get_value(Fields.Analysis.oetraj_rec)

            for fd in oetraj_rec.get_fields():

                name = fd.get_name()

                if name in "ProtTraj_OPLMD":
                    mdtrajrecord = MDDataRecord(oetraj_rec)
                    prot = mdtrajrecord.get_protein_traj
                    prot.SetTitle("ProteinTraj")
                    fn = os.path.join(title, "protein_trajectory.oeb")
                    with oechem.oemolostream(fn) as ofs:
                        oechem.OEWriteConstMolecule(ofs, prot)

                elif name in "WatTraj":
                    wat = oetraj_rec.get_value(fd)
                    wat.SetTitle("WatTraj")
                    fn = os.path.join(title, "water_trajectory.oeb")
                    with oechem.oemolostream(fn) as ofs:
                        oechem.OEWriteConstMolecule(ofs, wat)

                elif name in "LigTraj":
                    lig = oetraj_rec.get_value(fd)
                    lig.SetTitle("LigTraj")
                    fn = os.path.join(title, "ligand_trajectory.oeb")
                    with oechem.oemolostream(fn) as ofs:
                        oechem.OEWriteConstMolecule(ofs, lig)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):

        try:
            # Floe Report
            report_zip_fn = None
            if in_orion():
                session = OrionSession(
                    requests_session=get_session(
                        retry_dict={
                            403: 5,
                            404: 20,
                            409: 45,
                            460: 15,
                            500: 2,
                            502: 45,
                            503: 45,
                            504: 45,
                        }
                    )
                )

                report_zip_fn = "floe_report.zip"
                download_report(self.floe_report_collection_id, session, report_zip_fn)

            # Tar the files dir with its content:
            tar_fn = self.protein_name + '.tar.gz'

            with tarfile.open(tar_fn, mode='w:gz') as archive:
                for d in self.dirs:
                    archive.add(d)
                if report_zip_fn is not None:
                    archive.add(report_zip_fn)

            if in_orion():

                session = OrionSession(
                    requests_session=get_session(
                        retry_dict={
                            403: 5,
                            404: 20,
                            409: 45,
                            460: 15,
                            500: 2,
                            502: 45,
                            503: 45,
                            504: 45,
                        }
                    )
                )

                file_upload = File.upload(session, self.opt['out_file_name'], tar_fn)

                session.tag_resource(file_upload, "MD Data")

                job_id = environ.get('ORION_JOB_ID')

                if job_id:
                    session.tag_resource(file_upload, "Job {}".format(job_id))

            for d in self.dirs:
                shutil.rmtree(d)

            return

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.log.error(traceback.format_exc())


class ParallelClusterOETrajCube(ParallelMixin, ClusterOETrajCube):
    title = "Parallel " + ClusterOETrajCube.title
    description = "(Parallel) " + ClusterOETrajCube.description
    uuid = "216973c9-5f13-46f9-b79d-dee9d90398e9"


class ParallelMakeClusterTrajOEMols(ParallelMixin, MakeClusterTrajOEMols):
    title = "Parallel " + MakeClusterTrajOEMols.title
    description = "(Parallel) " + MakeClusterTrajOEMols.description
    uuid = "3415b9cb-ee9b-4138-a558-27e7f357d71d"


class ParallelClusterPopAnalysis(ParallelMixin, ClusterPopAnalysis):
    title = "Parallel " + ClusterPopAnalysis.title
    description = "(Parallel) " + ClusterPopAnalysis.description
    uuid = "bd607413-9fd8-4cb5-a802-54575cdf9833"


class ParallelTrajAnalysisReportDataset(ParallelMixin, TrajAnalysisReportDataset):
    title = "Parallel " + TrajAnalysisReportDataset.title
    description = "(Parallel) " + TrajAnalysisReportDataset.description
    uuid = "8d5308e6-58b4-4562-ae85-35bc47dcea0c"


class ParallelMDTrajAnalysisClusterReport(ParallelMixin,  MDTrajAnalysisClusterReport):
    title = "Parallel " + MDTrajAnalysisClusterReport.title
    description = "(Parallel) " + MDTrajAnalysisClusterReport.description
    uuid = "10f572c8-a874-47de-8f48-19ac76f72bdd"


