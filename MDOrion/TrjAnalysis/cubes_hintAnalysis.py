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

from MDOrion.Standards import Fields

import numpy as np
import json

import MDOrion.TrjAnalysis.utils as utl
import oetrajanalysis.trajOEHint_utils as hint

from openeye import oechem

import traceback

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from MDOrion.Standards.mdrecord import MDDataRecord


class BintScoreInitialPoseAndTrajectory(RecordPortsMixin, ComputeCube):
    title = 'Compute BintScore for Initial Pose and Trajectory'
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Binding Interactions', 'Ligand', 'Protein']

    description = """
    Compute BintScore for Initial Pose and Trajectory

    This Cube computes the BintScore for the Initial Pose and then 
    weights each of the BintScore interactions of the initial pose
    by their fractional occupancy over the selected frames of the
    ligand trajectory.
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
            opt['Logger'].info(' Beginning BintScoreInitialPoseAndTrajectory')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Compute BintScore for Initial Pose and Trajectory.'
                .format(system_title) )

            # Get the ligand which will be a multiconformer molecule with the starting
            # conformer for each simulation
            if not record.has_field(Fields.ligand):
                raise ValueError('{} could not find the ligand field'.format(system_title))
            ligand = utl.RequestOEFieldType(record, Fields.ligand)
            lig_name = utl.RequestOEFieldType(record, Fields.ligand_name)
            protein = utl.RequestOEFieldType(record, Fields.protein)

            # Get the ligand trajectory OEMol with one conformer per trajectory frame
            if not record.has_field(Fields.Analysis.oetraj_rec):
                raise ValueError('{} could not find the traj record'.format(system_title))
            opt['Logger'].info('{} found the traj record'.format(system_title))
            oetrajRecord = record.get_value(Fields.Analysis.oetraj_rec)
            ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} : got ligTraj with {} atoms, {} confs'.format(
                system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )

            # Extract the protTraj OEMol from the OETraj record
            mdtrajrecord = MDDataRecord(oetrajRecord)
            protTraj = mdtrajrecord.get_protein_traj
            opt['Logger'].info('{} got protTraj with {} atoms, {} confs'.format(
                system_title, protTraj.NumAtoms(), protTraj.NumConfs()) )


            # Generate HintDict of good hints and associated BintScore for initial pose
            good_hints_init_pose = hint.GoodHintDict(ligand, protein)

            # Calc BintScore for the initial pose
            initBintScore = hint.BintScoreFromHintDict(good_hints_init_pose,hint.intnStrengths)
            opt['Logger'].info('{} Initial pose for ligand {:s} has BintScore {:.1f}'.format(
                system_title, ligand.GetTitle(), initBintScore))

            # Calc list of per-frame BintScores for a trajectory
            trajBintScoreList = hint.TrajBintScoreListFromRefHints(ligTraj,protTraj,good_hints_init_pose)
            opt['Logger'].info('{} Traj has {} BintScores'.format(system_title,len(trajBintScoreList)))

            # Calc Trajectory BintScore: Initial BintScore weighted by fractional occupancy of each interaction
            trajBintScore, trajBintStderr, CI05, CI95 = utl.MeanAndBootstrapStdErrCI(trajBintScoreList)
            opt['Logger'].info('{} Trajectory of ligand {:s} has BintScore {:.2f} +/- {:.2f}'.format(
                system_title, ligTraj.GetTitle(), trajBintScore, trajBintStderr))

            # Create new record with Bint-related results
            bintRecord = OERecord()

            # Output results on the Bint record
            bintRecord.set_value(Fields.Bint.initBintScore, initBintScore)
            bintRecord.set_value(Fields.Bint.trajBintScoreList, trajBintScoreList)
            bintRecord.set_value(Fields.Bint.trajBintScore, trajBintScore)
            bintRecord.set_value(Fields.Bint.trajBintStderr, trajBintStderr)

            # The Bint record goes on the top-level record
            record.set_value(Fields.Bint.oebint_rec, bintRecord)

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ComparePoseBintsToTrajBintsCube on {}'.
                               format(str(e), system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class ParallelBintScoreInitialPoseAndTrajectory(ParallelMixin,  BintScoreInitialPoseAndTrajectory):
    title = "Parallel " + BintScoreInitialPoseAndTrajectory.title
    description = "(Parallel) " + BintScoreInitialPoseAndTrajectory.description
    #uuid = "10f572c8-a874-47de-8f48-19ac76f72bdd"

