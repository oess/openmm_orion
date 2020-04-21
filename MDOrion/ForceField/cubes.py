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

import traceback

from orionplatform.mixins import RecordPortsMixin

from floe.api import (ParallelMixin,
                      parameters,
                      ComputeCube)

from oeommtools import utils as oeommutils

from MDOrion.ForceField.ff_library import ff_library
from MDOrion.ForceField import utils

import parmed

from openeye import oechem

from oeommtools import data_utils as pack_utils

from MDOrion.Standards import MDStageNames, MDStageTypes

from MDOrion.Standards.mdrecord import MDDataRecord, Fields

from MDOrion.MDEngines.utils import MDState


from simtk.openmm import app

from simtk import unit

import os


class ForceFieldCube(RecordPortsMixin, ComputeCube):
    title = "Force Field Application"
    # version = "0.1.4"
    classification = [["Force Field"]]
    tags = ['ForceField']
    description = """
    This cube parametrizes a flask with the selected force fields. 
    The cube tries to split a flask into components: protein, ligand, 
    water and excipients. The user can select the parametrization to be 
    applied to each component. The protein forcefield is limited to 
    standard amino acids and limited support to non-standard. Sugars 
    are not currently supported but this will be improved in coming 
    releases. The cube requires a record as input and produces a new 
    record where the flask has been parametrized. The parametrization 
    is carried out by using a Parmed object 
    (https://github.com/ParmEd/ParmEd) 
    which will be present on the emitted record. The supported protein 
    force fields are amber99sb-ildn and the new amberfb-15. Small organic
    molecules like ligands and excipients can be parametrized by using 
    GAFF, GAFF2 and SMIRNOFF forcefields. The flask splitting is based on the ligand 
    residue name. The default one is “LIG” and can be changed by using 
    the provided cube parameter. Water is currently parametrized by 
    using TIP3P force field water model only.
    """

    uuid = "aac0d06f-afd3-4801-ba50-2d703a07ab35"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_forcefield = parameters.StringParameter(
        'protein_forcefield',
        default=sorted(ff_library.proteinff)[0],
        choices=sorted(ff_library.proteinff),
        help_text='Force field parameters to be applied to the protein')

    solvent_forcefield = parameters.StringParameter(
        'solvent_forcefield',
        default=sorted(ff_library.solventff)[0],
        help_text='Force field parameters to be applied to the water')

    ligand_forcefield = parameters.StringParameter(
        'ligand_forcefield',
        default=sorted(ff_library.ligandff)[0],
        choices=sorted(ff_library.ligandff),
        help_text='Force field to be applied to the ligand')

    suffix = parameters.StringParameter(
        'suffix',
        default='prep',
        help_text='Filename suffix for output simulation files')

    other_forcefield = parameters.StringParameter(
        'other_forcefield',
        default=sorted(ff_library.otherff)[0],
        choices=sorted(ff_library.otherff),
        help_text='Force field used to parametrize other molecules not recognized by the '
                  'protein force field like excipients')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            opt['CubeTitle'] = self.title

            if not record.has_value(Fields.design_unit):
                raise ValueError("The Design Unit is not present on the record")

            du = record.get_value(Fields.design_unit)

            # Design Unit Parametrization
            utils.parametrize_du(du, opt)











            # # Create the MD record to use the MD Record API
            # mdrecord = MDDataRecord(record)
            #
            # flask = mdrecord.get_flask
            #
            # if not mdrecord.has_title:
            #     self.log.warn("Missing record Title field")
            #     flask_title = flask.GetTitle()[0:12]
            # else:
            #     flask_title = mdrecord.get_title
            #
            # # Split the complex in components in order to apply the FF
            # protein, ligand, water, excipients = oeommutils.split(flask, ligand_res_name=opt['lig_res_name'])
            #
            # self.log.info("[{}] Components of flask {}:\n  Protein atoms = {}\n  Ligand atoms = {}\n"
            #               "  Water atoms = {}\n  Excipients atoms = {}".format(opt['CubeTitle'],
            #                                                                    flask_title,
            #                                                                    protein.NumAtoms(),
            #                                                                    ligand.NumAtoms(),
            #                                                                    water.NumAtoms(),
            #                                                                    excipients.NumAtoms()))
            # protein = ffutils.clean_tags(protein)
            # ligand = ffutils.clean_tags(ligand)
            # water = ffutils.clean_tags(water)
            # excipients = ffutils.clean_tags(excipients)
            #
            # sys_id = mdrecord.get_flask_id
            #
            # # Unique prefix name used to output parametrization files
            # opt['prefix_name'] = flask_title + '_'+str(sys_id)
            #
            # oe_mol_list = []
            # par_mol_list = []
            #
            # # Apply FF to the Protein
            # if protein.NumAtoms():
            #     oe_mol_list.append(protein)
            #     protein_structure = ffutils.applyffProtein(protein, opt)
            #     par_mol_list.append(protein_structure)
            #
            # # Apply FF to the ligand
            # if ligand.NumAtoms():
            #     oe_mol_list.append(ligand)
            #     ligand_structure = ffutils.applyffLigand(ligand, opt)
            #     par_mol_list.append(ligand_structure)
            #
            # # Apply FF to water molecules
            # if water.NumAtoms():
            #     oe_mol_list.append(water)
            #     water_structure = ffutils.applyffWater(water, opt)
            #     par_mol_list.append(water_structure)
            #
            # # Apply FF to the excipients
            # if excipients.NumAtoms():
            #     oe_mol_list.append(excipients)
            #     excipient_structure = ffutils.applyffExcipients(excipients, opt)
            #     par_mol_list.append(excipient_structure)
            #
            # # Build the overall Parmed structure
            # flask_structure = parmed.Structure()
            #
            # for struc in par_mol_list:
            #     flask_structure = flask_structure + struc
            #
            # flask_reassembled = oe_mol_list[0].CreateCopy()
            # num_atom_flask = flask_reassembled.NumAtoms()
            #
            # for idx in range(1, len(oe_mol_list)):
            #     oechem.OEAddMols(flask_reassembled, oe_mol_list[idx])
            #     num_atom_flask += oe_mol_list[idx].NumAtoms()
            #
            # if not num_atom_flask == flask_structure.topology.getNumAtoms():
            #     raise ValueError("Parmed and OE topologies mismatch atom number {} vs {}"
            #                      .format(num_atom_flask, flask_structure.topology.getNumAtoms()))
            #
            # flask_reassembled.SetTitle(flask.GetTitle())
            #
            # # Set Parmed structure box_vectors
            # is_periodic = True
            #
            # try:
            #     vec_data = pack_utils.getData(flask_reassembled, tag='box_vectors')
            #     vec = pack_utils.decodePyObj(vec_data)
            #     flask_structure.box_vectors = vec
            # except:
            #     is_periodic = False
            #     self.log.warn("flask {} has been parametrize without periodic box vectors ".format(flask_title))
            #
            # # Set atom serial numbers, Ligand name and HETATM flag
            # for at in flask_reassembled.GetAtoms():
            #     thisRes = oechem.OEAtomGetResidue(at)
            #     thisRes.SetSerialNumber(at.GetIdx())
            #     if thisRes.GetName() == 'UNL':
            #         # thisRes.SetName("LIG")
            #         thisRes.SetHetAtom(True)
            #     oechem.OEAtomSetResidue(at, thisRes)
            #
            # if flask_reassembled.NumAtoms() != flask_structure.topology.getNumAtoms():
            #     raise ValueError("OEMol flask {} and generated Parmed structure "
            #                      "mismatch atom numbers".format(flask_title))
            #
            # flask_formal_charge = 0
            # for at in flask_reassembled.GetAtoms():
            #     flask_formal_charge += at.GetFormalCharge()
            #
            # flask_partial_charge = 0.0
            # for at in flask_structure.atoms:
            #     flask_partial_charge += at.charge
            #
            # if abs(flask_formal_charge - flask_partial_charge) > 0.01:
            #     raise ValueError("flask Formal charge and flask Partial charge mismatch: {} vs {}".format(
            #         flask_formal_charge, flask_partial_charge))
            #
            # # Copying the charges between the parmed structure and the oemol
            # for parm_at, oe_at in zip(flask_structure.atoms, flask_reassembled.GetAtoms()):
            #
            #     if parm_at.atomic_number != oe_at.GetAtomicNum():
            #         raise ValueError("Atomic number mismatch between the Parmed and the OpenEye topologies: {} - {}".
            #                          format(parm_at.atomic_number, oe_at.GetAtomicNum()))
            #
            #     oe_at.SetPartialCharge(parm_at.charge)
            #
            # # Check if it is possible to create the OpenMM System
            # if is_periodic:
            #     omm_flask = flask_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
            #                                              nonbondedCutoff=10.0 * unit.angstroms,
            #                                              constraints=None,
            #                                              removeCMMotion=False,
            #                                              rigidWater=False)
            # else:
            #     omm_flask = flask_structure.createSystem(nonbondedMethod=app.NoCutoff,
            #                                              constraints=None,
            #                                              removeCMMotion=False,
            #                                              rigidWater=False)
            # mdrecord.set_title(flask_title)
            # mdrecord.set_flask(flask_reassembled)
            #
            # mdrecord.set_parmed(flask_structure, shard_name="Parmed_" + flask_title + '_' + str(sys_id))
            #
            # data_fn = os.path.basename(mdrecord.cwd) + '_' + flask_title+'_' + str(sys_id) + '-' + opt['suffix']+'.tar.gz'
            #
            # if not mdrecord.add_new_stage(MDStageNames.ForceField,
            #                               MDStageTypes.SETUP,
            #                               flask_reassembled,
            #                               MDState(flask_structure),
            #                               data_fn):
            #     raise ValueError("Problems adding the new Parametrization Stage")
            #
            # self.success.emit(mdrecord.get_record)
            #
            # del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class ParallelForceFieldCube(ParallelMixin, ForceFieldCube):
    title = "Parallel " + ForceFieldCube.title
    description = "(Parallel) " + ForceFieldCube.description
    uuid = "deb6b453-0ddf-4f1c-a709-cda1f3c47af1"

