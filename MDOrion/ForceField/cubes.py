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
from floe.api import (ParallelMixin,
                      parameter)

from cuberecord import OERecordComputeCube

from oeommtools import utils as oeommutils

from MDOrion.ForceField import utils as ffutils

import parmed

from openeye import oechem

from oeommtools import data_utils as pack_utils

from MDOrion.Standards import MDStageTypes

from MDOrion.Standards.mdrecord import MDDataRecord

from MDOrion.MDEngines.utils import MDState

from simtk.openmm import app

from simtk import unit

import os


class ForceFieldCube(ParallelMixin, OERecordComputeCube):
    title = "Force Field Application"
    version = "0.1.0"
    classification = [["Force Field"]]
    tags = ['ForceField']
    description = """
    This cube parametrized a system with the selected force fields. 
    The cube tries to split a system into components: protein, ligand, 
    water and excipients. The user can select the parametrization to be 
    applied to each component. The protein forcefield is limited to 
    standard amino acids and limited support to non-standard. Sugars 
    are not currently supported but this will be improved in coming 
    releases. The cube requires a record as input and produces a new 
    record where the system has been parametrized. The parametrization 
    is carried out by using a Parmed object 
    (https://github.com/ParmEd/ParmEd) 
    which will be present on the emitted record. The supported protein 
    force fields are amber99sb-ildn and the new amberfb-15. Small organic
    molecules like ligands and excipients can be parametrized by using 
    GAFF, GAFF2 and SMIRNOFF. The system spitting is based on the ligand 
    residue name. The default one is “LIG” and can be changed by using 
    the provided cube parameter. Water is currently parametrized by 
    using TIP3P force field water model only.

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the systems to parametrize

    Output:
    -------
    oechem.OEDataRecord - Streamed-out of records with the parametrized systems.
    Each record will contain a new Parmed object that carry out the 
    system parametrization
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_forcefield = parameter.StringParameter(
        'protein_forcefield',
        default=sorted(ffutils.proteinff)[0],
        choices=sorted(ffutils.proteinff),
        help_text='Force field parameters to be applied to the protein')

    solvent_forcefield = parameter.StringParameter(
        'solvent_forcefield',
        default=sorted(ffutils.solventff)[0],
        help_text='Force field parameters to be applied to the water')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        default=sorted(ffutils.ligandff)[0],
        choices=sorted(ffutils.ligandff),
        help_text='Force field to be applied to the ligand')

    lig_res_name = parameter.StringParameter(
        'lig_res_name',
        default='LIG',
        help_text='Ligand residue name. This is used during the spitting to identify the ligand')

    suffix = parameter.StringParameter(
        'suffix',
        default='prep',
        help_text='Filename suffix for output simulation files')

    other_forcefield = parameter.StringParameter(
        'other_forcefield',
        default=sorted(ffutils.otherff)[0],
        choices=sorted(ffutils.otherff),
        help_text='Force field used to parametrize other molecules not recognized by the '
                  'protein force field like excipients')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            opt['CubeTitle'] = self.title

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            system = mdrecord.get_primary

            if not mdrecord.has_title:
                self.log.warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = mdrecord.get_title

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = oeommutils.split(system, ligand_res_name=opt['lig_res_name'])

            self.log.info("[{}] \nComplex name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(opt['CubeTitle'],
                                                                                         system_title,
                                                                                         protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))
            sys_id = mdrecord.get_id

            # Unique prefix name used to output parametrization files
            opt['prefix_name'] = system_title + '_'+str(sys_id)

            oe_mol_list = []
            par_mol_list = []

            # Apply FF to the Protein
            if protein.NumAtoms():
                oe_mol_list.append(protein)
                protein_structure = ffutils.applyffProtein(protein, opt)
                par_mol_list.append(protein_structure)

            # Apply FF to the ligand
            if ligand.NumAtoms():
                oe_mol_list.append(ligand)
                ligand_structure = ffutils.applyffLigand(ligand, opt)
                par_mol_list.append(ligand_structure)

            # Apply FF to water molecules
            if water.NumAtoms():
                oe_mol_list.append(water)
                water_structure = ffutils.applyffWater(water, opt)
                par_mol_list.append(water_structure)

            # Apply FF to the excipients
            if excipients.NumAtoms():
                excipient_structure = ffutils.applyffExcipients(excipients, opt)
                par_mol_list.append(excipient_structure)

                # The excipient order is set equal to the order in related
                # parmed structure to avoid possible atom index mismatching
                excipients = oeommutils.openmmTop_to_oemol(excipient_structure.topology,
                                                           excipient_structure.positions,
                                                           verbose=False)
                oechem.OEPerceiveBondOrders(excipients)
                oe_mol_list.append(excipients)

            # Build the overall Parmed structure
            system_structure = parmed.Structure()

            for struc in par_mol_list:
                system_structure = system_structure + struc

            system_reassembled = oe_mol_list[0].CreateCopy()
            num_atom_system = system_reassembled.NumAtoms()

            for idx in range(1, len(oe_mol_list)):
                oechem.OEAddMols(system_reassembled, oe_mol_list[idx])
                num_atom_system += oe_mol_list[idx].NumAtoms()

            if not num_atom_system == system_structure.topology.getNumAtoms():
                raise ValueError("Parmed and OE topologies mismatch atom number "
                                 "error for system: {}".format(system_title))

            system_reassembled.SetTitle(system.GetTitle())

            # Set Parmed structure box_vectors
            is_periodic = True

            try:
                vec_data = pack_utils.getData(system_reassembled, tag='box_vectors')
                vec = pack_utils.decodePyObj(vec_data)
                system_structure.box_vectors = vec
            except:
                is_periodic = False
                self.log.warn("System {} has been parametrize without periodic box vectors ".format(system_title))

            # Set atom serial numbers, Ligand name and HETATM flag
            for at in system_reassembled.GetAtoms():
                thisRes = oechem.OEAtomGetResidue(at)
                thisRes.SetSerialNumber(at.GetIdx())
                if thisRes.GetName() == 'UNL':
                    # thisRes.SetName("LIG")
                    thisRes.SetHetAtom(True)
                oechem.OEAtomSetResidue(at, thisRes)

            if system_reassembled.GetMaxAtomIdx() != system_structure.topology.getNumAtoms():
                raise ValueError("OEMol system {} and generated Parmed structure "
                                 "mismatch atom numbers".format(system_title))

            # Copying the charges between the parmed structure and the oemol
            for parm_at, oe_at in zip(system_structure.atoms, system_reassembled.GetAtoms()):

                if parm_at.atomic_number != oe_at.GetAtomicNum():
                    raise ValueError("Atomic number mismatch between the Parmed and the OpenEye topologies: {} - {}".
                                     format(parm_at.atomic_number, oe_at.GetAtomicNum()))

                oe_at.SetPartialCharge(parm_at.charge)

            # Check if it is possible to create the OpenMM System
            if is_periodic:
                omm_system = system_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                                           nonbondedCutoff=10.0 * unit.angstroms,
                                                           constraints=None,
                                                           removeCMMotion=False,
                                                           rigidWater=False)
            else:
                omm_system = system_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                           constraints=None,
                                                           removeCMMotion=False,
                                                           rigidWater=False)
            mdrecord.set_title(system_title)
            mdrecord.set_primary(system_reassembled)

            mdrecord.set_parmed(system_structure, shard_name="Parmed_" + system_title + '_' + str(sys_id))

            data_fn = os.path.basename(mdrecord.cwd) + '_' + system_title+'_' + str(sys_id) + '-' + opt['suffix']+'.tar.gz'

            if not mdrecord.add_new_stage(self.title,
                                          MDStageTypes.SETUP,
                                          system_reassembled,
                                          MDState(system_structure),
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
