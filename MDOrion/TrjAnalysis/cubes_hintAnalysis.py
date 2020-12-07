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

import MDOrion.TrjAnalysis.utils as utl
import MDOrion.TrjAnalysis.trajOEHint_utils as hint

from openeye import oechem

import traceback

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from MDOrion.Standards.mdrecord import MDDataRecord


class ComparePoseBintsToTrajBints(RecordPortsMixin, ComputeCube):
    title = 'Compare Binding Interations of Pose To Traj'
    # version = "0.1.4"
    classification = [["Analysis"]]
    tags = ['Binding Interactions', 'Ligand', 'Protein']

    description = """
    Compare Binding Interations of Pose To Traj

    This Cube compares the initial pose of the ligand to the
    ligand trajectory in terms of the OEHint binding interactions.
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
            opt['Logger'].info(' Beginning ComparePoseBintsToTrajBintsCube')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to compare initial and traj Bints.'
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


            # Binding interactions (OEHints) for the initial pose

            # Perceive OEHints
            asite = oechem.OEInteractionHintContainer(protein, ligand)
            if not oechem.OEIsValidActiveSite(asite):
                raise ValueError("{} Cannot initialize OEInteractionHintContainer".format(system_title))
            oechem.OEPerceiveInteractionHints(asite)

            # Make a dict of good interactions and calculate initBintScore
            good_hints_init_pose = hint.DictAllowedInteractionStrings(asite, hint.goodIntnTypes)
            initBintScore = 0
            for interType in good_hints_init_pose.keys():
                numInters = len(good_hints_init_pose[interType])
                #print(interType, numInters)
                initBintScore -= float(numInters)
            opt['Logger'].info('{} Initial pose for ligand {:s} gets a Bint score of {:.2f}'.
                               format(system_title,ligand.GetTitle(), initBintScore))

            # Make dict-of-dicts for holding count of traj occurrence for each interaction
            hint_traj = dict()
            for key in good_hints_init_pose.keys():
                #print(key)
                inter_dict = dict()
                for inter in good_hints_init_pose[key]:
                    inter_dict[inter] = 0
                hint_traj[key] = inter_dict


            # Binding interactions (OEHints) for each fram of the trajectory

            # Loop over lig, prot confs (traj frames)), calculating HintInteractions and
            # counting overlap of favorable hints with the initial pose
            for cLig, cProt in zip(ligTraj.GetConfs(), protTraj.GetConfs()):
                frameid = cLig.GetIdx()

                # Perceive interactions on paired conformers of protein and ligand
                frameSite = oechem.OEInteractionHintContainer(oechem.OEMol(cProt), oechem.OEMol(cLig))
                if not oechem.OEIsValidActiveSite(frameSite):
                    raise ValueError("{} Cannot initialize active site of paired conformers of protein and ligand!".
                                     format(system_title))
                oechem.OEPerceiveInteractionHints(frameSite)

                # Make a dict of good interactions
                good_hints_frame = hint.DictAllowedInteractionStrings(frameSite, hint.goodIntnTypes)

                # find conformer good hints that match init_pose good hints
                if not hint.AddFrameHintOverlapToRefHints(hint_traj, good_hints_frame, frameid):
                    raise ValueError("{} Unsuccessful with frame {}".format(system_title,frameid))

            # Calculate trajBintScore
            nConfs = ligTraj.NumConfs()
            trajBintScore = 0.0
            for interType in hint_traj.keys():
                #print(interType)
                for inter in hint_traj[interType].keys():
                    fracOcc = float(hint_traj[interType][inter] / nConfs)
                    trajBintScore -= fracOcc
                    #print('  <{:s}> : {:.3f}'.format(inter, fracOcc))
            trajBintOcc = trajBintScore/initBintScore
            opt['Logger'].info('{} ligand {:s} traj maintains {:.4f} of initial good binding interactions'.
                  format(system_title,ligTraj.GetTitle(),trajBintOcc))
            opt['Logger'].info('{} ligand {:s} trajectory gets a Bint score of {:.2f}'.
                  format(system_title,ligTraj.GetTitle(), trajBintScore))


            # Output results on the top-level record
            record.set_value(OEField('InitBintScore', Types.Float), initBintScore)
            record.set_value(OEField('TrajBintOcc', Types.Float), trajBintOcc)
            record.set_value(OEField('TrajBintScore', Types.Float), trajBintScore)

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in ComparePoseBintsToTrajBintsCube on {}'.
                               format(str(e), system_title))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class ParallelComparePoseBintsToTrajBints(ParallelMixin,  ComparePoseBintsToTrajBints):
    title = "Parallel " + ComparePoseBintsToTrajBints.title
    description = "(Parallel) " + ComparePoseBintsToTrajBints.description
    #uuid = "10f572c8-a874-47de-8f48-19ac76f72bdd"

