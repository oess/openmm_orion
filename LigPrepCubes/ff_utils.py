import subprocess, tempfile, parmed
from openeye import oechem, oequacpac
import openmoltools
from openmoltools.openeye import *


from big_storage import LargeFile
from floe.api.orion import in_orion


def assignELF10charges(molecule, max_confs=800, strictStereo=True):
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

    # ELF10 charges
    status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCELF10Charges())

    if not status:
        raise RuntimeError("OEAssignCharges returned error code %d" % status)

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
        if not forcefield in ['SMIRNOFF', 'GAFF', 'GAFF2']:
            raise RuntimeError('Selected forcefield %s is not GAFF/GAFF2/SMIRNOFF' % forcefield)
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
            cmd.write( "source leaprc.%s; quit" % self.forcefield.lower() )
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
        from openforcefield.utils.utils import generateSMIRNOFFStructure
        if not molecule:
            molecule = self.molecule
        try:
            molecule_structure = generateSMIRNOFFStructure(molecule)
        except:
            raise RuntimeError('Error generating SMIRNOFF Structure for %s' % molecule.GetTitle())
        return molecule_structure

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

        # Load via ParmEd
        molecule_structure = parmed.amber.AmberParm(prmtop, inpcrd)

        if self.delete_out_files:
            os.remove(gaff_mol2_filename)
            os.remove(frcmod_filename)
            os.remove(prmtop)
            os.remove(inpcrd)

        return molecule_structure

    def parameterize(self):
        if self.forcefield == 'SMIRNOFF':
            structure = self.getSmirnoffStructure()
        elif self.forcefield in ['GAFF', 'GAFF2']:
            structure = self.getGaffStructure()
        self.structure = structure
        return self.structure


def upload(filename):

    file_id = filename

    if in_orion():
        file_id = LargeFile(filename)

    return file_id


def download(file_id):

    filename = file_id

    if in_orion():
        filename = file_id.retrieve()
        file_id.delete()

    return filename