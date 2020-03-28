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


import subprocess
import tempfile
import parmed

from openeye import oechem, oequacpac
import openmoltools
from openmoltools.openeye import *

from simtk.openmm.app import AmberInpcrdFile, AmberPrmtopFile
from simtk.openmm import app

from openforcefield.typing.engines.smirnoff import ForceField

from pkg_resources import resource_filename

from openforcefield.topology import Topology, Molecule

from oeommtools.utils import oemol_to_openmmTop

ligandff = {'Gaff': 'GAFF',
            'Gaff2': 'GAFF2',
            'Smirnoff99Frosst': 'smirnoff99Frosst.offxml',
            'OpenFF_1.1.0': "openff_unconstrained-1.1.0.offxml"}


def assignELF10charges(molecule, max_confs=800, strictStereo=True, opt=None):
    """
     This function computes atomic partial charges for an OEMol by
     using the ELF10 method

    Parameters:
    -----------
    molecule : OEMol object
        The molecule that needs to be charged
    max_confs : integer
        The max number of conformers used to calculate the atomic partial charges
    strictStereo : bool
        a flag used to check if atoms need to have assigned stereo chemistry or not

    Return:
    -------
    mol_copy : OEMol
        a copy of the original molecule with assigned atomic partial charges
    """

    mol_copy = molecule.CreateCopy()

    # The passed molecule could have already conformers. If the conformer number
    # does not exceed the max_conf threshold then max_confs conformations will
    # be generated
    if not mol_copy.GetMaxConfIdx() > 200:
        # Generate up to max_confs conformers
        mol_copy = generate_conformers(mol_copy, max_confs=max_confs, strictStereo=strictStereo)

    # Assign MMFF Atom types
    if not oechem.OEMMFFAtomTypes(mol_copy):
        raise RuntimeError("MMFF atom type assignment returned errors")

    # Check for Carboxylic Acid patterns in the molecule
    smarts = '(O=)[C][O,S][H]'
    ss = oechem.OESubSearch(smarts)

    oechem.OEPrepareSearch(mol_copy, ss)
    unique_match = True

    a_match_list = []
    for match in ss.Match(mol_copy, unique_match):

        for ma in match.GetAtoms():
            a_match_list.append(ma.target)

    # Set the Carboxylic Acid torsion to zero for each generated conformers
    if a_match_list:

        if len(a_match_list) % 4 != 0:
            raise ValueError("The atom matching list must be multiple of 4")

        for i in range(0, len(a_match_list), 4):

            chunk = a_match_list[i:i + 4]

            for conf in mol_copy.GetConfs():

                conf.SetTorsion(chunk[0],
                                chunk[1],
                                chunk[2],
                                chunk[3], 0.0)

    # Try to calculate the ELF10 charges for the molecule
    quacpac_status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCELF10Charges())

    if not quacpac_status:
        opt['Logger'].warn("OEAM1BCCELF10 charge assignment failed downgrading "
                           "to OEAM1BCC charge assignment for this molecule: {}".format(mol_copy.GetTitle()))

        quacpac_status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCCharges())

    if not quacpac_status:
        raise RuntimeError("OEAssignCharges returned error code {}".format(quacpac_status))

    return mol_copy


class ParamLigStructure(object):
    """
    Generates parametrized ParmEd structure of the molecule with a chosen force field

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The openeye molecule to be parameterized
    forcefield : str
        String specifying the forcefield parameters to be used
    prefix_name : str
        String specifying the output prefix filename

    Returns
    ---------
    packedmol : openeye.oechem.OEMol
        Openeye molecule with the ParmEd Structure attached.
    """

    def __init__(self, molecule, forcefield, prefix_name='ligand', delete_out_files=True):
        if not forcefield in list(ligandff.values()):
            raise RuntimeError('The selected ligand force field is not '
                               'supported {}. Available {}'.format(forcefield, list(ligandff.keys())))
        else:
            self.molecule = molecule
            self.forcefield = str(forcefield).strip()
            self.structure = None
            self.prefix_name = prefix_name
            self.delete_out_files = delete_out_files

    @staticmethod
    def checkTleap(self):
        # Try to check if tleap is going to fail
        with open('tleap_commands', 'w') as cmd:
            cmd.write("source leaprc.%s; quit" % self.forcefield.lower())
        tmp = subprocess.getoutput('tleap -f tleap_commands')
        elements = tmp.split('\n')
        for elem in elements:
            if 'Could not open file' in elem:
                raise RuntimeError('Error encountered trying to load %s in tleap.'% self.forcefield)
            if 'command not found' in elem:
                raise RuntimeError('Error: requires tleap.')
        return True

    def checkCharges(self, molecule):
        # Check that molecule is charged.
        is_charged = False
        for atom in molecule.GetAtoms():
            if atom.GetPartialCharge() != 0.0:
                is_charged = True
        if not is_charged:
            raise Exception('Molecule %s has no charges; input molecules must be charged.' % molecule.GetTitle())

    def getSmirnoffStructure(self, molecule=None):

        if not molecule:
            molecule = self.molecule
        try:
            self.checkCharges(molecule)
        except:
            print("WARNING: Missing Charges, assigning elf10 charges to molecule")

            molecule = assignELF10charges(molecule)

        if self.forcefield == ligandff['Smirnoff99Frosst']:

            fffn = resource_filename('openforcefield', os.path.join('data', 'test_forcefields/' + self.forcefield))

            if not os.path.exists(fffn):
                raise ValueError(
                    "Sorry! {} does not exist. If you just added it, you'll have to re-install".format(fffn))

            with open(fffn) as ffxml:
                ff = ForceField(ffxml, allow_cosmetic_attributes=True)

        elif self.forcefield == ligandff['OpenFF_1.1.0']:

            ff = ForceField(self.forcefield, allow_cosmetic_attributes=True)

        else:
            raise ValueError("Force Field not Supported: {}".format(self.forcefield))

        mol_off = Molecule.from_openeye(molecule, allow_undefined_stereo=True)
        topology = Topology.from_molecules([mol_off])

        omm_sys = ff.create_openmm_system(topology, charge_from_molecules=[mol_off])

        # omm_top = generateTopologyFromOEMol(molecule)
        # positions = mol_off.conformers[0]

        omm_top, positions = oemol_to_openmmTop(molecule)

        pmd_structure = parmed.openmm.load_topology(omm_top, omm_sys, xyz=positions)

        return pmd_structure

    def getGaffStructure(self, molecule=None, forcefield=None):
        if not molecule:
            molecule = self.molecule
        if not forcefield:
            forcefield = self.forcefield

        #  Try to check if tleap is going to fail
        self.checkCharges(molecule)

        # Determine formal charge (antechamber needs as argument)
        chg = 0
        for atom in molecule.GetAtoms():
            chg += atom.GetFormalCharge()

        # Write out mol to a mol2 file to process via AmberTools
        mol2file = tempfile.NamedTemporaryFile(suffix='.mol2')
        mol2filename = mol2file.name
        with oechem.oemolostream(mol2filename) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, molecule)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing molecule %s to mol2." % molecule.GetTitle())

        # Run antechamber to type and parmchk for frcmod
        # requires openmoltools 0.7.5 or later, which should be conda-installable via omnia
        gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber(self.prefix_name, mol2filename,
                                                                                 gaff_version=forcefield.lower(),
                                                                                 charge_method=None)

        # Run tleap using specified forcefield
        prmtop, inpcrd = openmoltools.amber.run_tleap(self.prefix_name, gaff_mol2_filename,
                                                      frcmod_filename,
                                                      leaprc='leaprc.%s' % forcefield.lower())

        # TODO Load via ParmEd: This is causing Problems
        #  Merging two structures (OpenMM PMD structure and
        #  Amber PMD Structure): The NB exception list is messed up
        # molecule_structure = parmed.amber.AmberParm(prmtop, inpcrd)

        # TODO MODIFIED BY GAC
        omm_prmtop = AmberPrmtopFile(prmtop)
        omm_inpcrd = AmberInpcrdFile(inpcrd)

        omm_system = omm_prmtop.createSystem(nonbondedMethod=app.NoCutoff)

        molecule_structure = parmed.openmm.load_topology(omm_prmtop.topology, omm_system, xyz=omm_inpcrd.positions)

        if self.delete_out_files:
            os.remove(gaff_mol2_filename)
            os.remove(frcmod_filename)
            os.remove(prmtop)
            os.remove(inpcrd)

        return molecule_structure

    def parameterize(self):

        if self.forcefield in [ligandff['OpenFF_1.1.0'], ligandff['Smirnoff99Frosst']]:
            structure = self.getSmirnoffStructure()

        elif self.forcefield in ['GAFF', 'GAFF2']:
            structure = self.getGaffStructure()

        self.structure = structure
        return self.structure

