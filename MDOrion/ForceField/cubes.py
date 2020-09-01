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

import traceback

from orionplatform.mixins import RecordPortsMixin

from floe.api import (ParallelMixin,
                      parameters,
                      ComputeCube)

from oemdtoolbox.ForceField.ff_library import ff_library

from MDOrion.Standards import MDStageNames, MDStageTypes

from MDOrion.Standards.mdrecord import MDDataRecord, Fields

from MDOrion.MDEngines.utils import MDState

from openeye import oechem

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

    ligand_forcefield = parameters.StringParameter(
        'ligand_forcefield',
        default=sorted(ff_library.ligandff)[0],
        choices=sorted(ff_library.ligandff),
        help_text='Force field to be applied to the ligand')

    suffix = parameters.StringParameter(
        'suffix',
        default='prep',
        help_text='Filename suffix for output simulation files')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            opt['CubeTitle'] = self.title

            if not record.has_value(Fields.md_components):
                raise ValueError("MD Components Field is missing")

            md_components = record.get_value(Fields.md_components)
            flask, map_comp = md_components.create_flask

            opt['Logger'].info(md_components.get_info)

            if not record.has_value(Fields.title):
                self.log.warn("Missing record Title field")
                flask_title = flask.GetTitle()[0:12]
            else:
                flask_title = record.get_value(Fields.title)

            # Parametrize the whole flask
            flask_pmd_structure = md_components.parametrize_components(protein_ff=opt['protein_forcefield'],
                                                                       ligand_ff=opt['ligand_forcefield'])

            # Set Parmed structure box_vectors
            is_periodic = True
            if md_components.get_box_vectors is not None:
                flask_pmd_structure.box_vectors = md_components.get_box_vectors
            else:
                is_periodic = False
                self.log.warn("Flask {} has been parametrize without periodic box vectors ".format(flask_title))

            if flask.NumAtoms() != flask_pmd_structure.topology.getNumAtoms():
                raise ValueError("The flask {} and the generated Parmed structure "
                                 "have mismatch atom numbers: {} vs {}".
                                 format(flask_title,
                                        flask.NumAtoms(),
                                        flask_pmd_structure.topology.getNumAtoms()))

            # Check Formal vs Partial charges
            flask_formal_charge = 0
            for at in flask.GetAtoms():
                flask_formal_charge += at.GetFormalCharge()

            flask_partial_charge = 0.0
            for at in flask_pmd_structure.atoms:
                flask_partial_charge += at.charge

            if abs(flask_formal_charge - flask_partial_charge) > 0.01:
                raise ValueError("Flask Formal charge and flask Partial charge mismatch: {} vs {}".format(
                    flask_formal_charge, flask_partial_charge))

            # Copying the charges between the parmed structure and the oemol
            for parm_at, oe_at in zip(flask_pmd_structure.atoms, flask.GetAtoms()):

                if parm_at.atomic_number != oe_at.GetAtomicNum():
                    raise ValueError(
                        "Atomic number mismatch between the Parmed and the OpenEye topologies: {} - {}".
                        format(parm_at.atomic_number, oe_at.GetAtomicNum()))

                oe_at.SetPartialCharge(parm_at.charge)

            # Set the component charges
            for comp_name, comp in md_components.get_components.items():
                for at_comp in comp.GetAtoms():

                    pred = oechem.OEHasAtomIdx(map_comp[comp_name][at_comp.GetIdx()])
                    at_flask = flask.GetAtom(pred)

                    if at_flask.GetAtomicNum() != at_comp.GetAtomicNum():
                        "Atomic number mismatch between the component {}  atom {} and the flask atom {}".\
                            format(comp_name, at_comp, at_flask)

                    at_comp.SetPartialCharge(at_flask.GetPartialCharge())

                md_components.set_component_by_name(comp_name, comp)

            # Update the components after setting the charges
            record.set_value(Fields.md_components, md_components)

            # Check if it is possible to create the OpenMM System
            if is_periodic:
                omm_flask = flask_pmd_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                                             nonbondedCutoff=10.0 * unit.angstroms,
                                                             constraints=None,
                                                             removeCMMotion=False,
                                                             rigidWater=False)
            else:
                omm_flask = flask_pmd_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                             constraints=None,
                                                             removeCMMotion=False,
                                                             rigidWater=False)

            mdrecord = MDDataRecord(record)
            sys_id = mdrecord.get_flask_id
            mdrecord.set_flask(flask)

            mdrecord.set_parmed(flask_pmd_structure, shard_name="Parmed_" + flask_title + '_' + str(sys_id))

            data_fn = os.path.basename(mdrecord.cwd) + '_' + flask_title+'_' + str(sys_id) + '-' + opt['suffix']+'.tar.gz'

            if not mdrecord.add_new_stage(MDStageNames.ForceField,
                                          MDStageTypes.SETUP,
                                          flask,
                                          MDState(flask_pmd_structure),
                                          data_fn):
                raise ValueError("Problems adding the new Parametrization Stage")

            self.success.emit(mdrecord.get_record)

            del mdrecord

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

