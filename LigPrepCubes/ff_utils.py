import io, os, traceback, string, random, subprocess, tempfile, parmed
from openeye import oechem, oequacpac
import openmoltools
from openmoltools.openeye import *
from OpenMMCubes.utils import get_data_filename


def assignCharges(molecule, max_confs=800, strictStereo=True, normalize=True, keep_confs=None):
    """Generate charges for an OpenEye OEMol molecule.
    Adopted get_charges() from
    https://github.com/choderalab/openmoltools/blob/master/openmoltools/openeye.py
    to use new oequacpac.OEAssignCharges()
    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers.
        Omega will be used to generate max_confs conformations.
    max_confs : int, optional, default=800
        Max number of conformers to generate
    strictStereo : bool, optional, default=True
        If False, permits smiles strings with unspecified stereochemistry.
        See https://docs.eyesopen.com/omega/usage.html
    normalize : bool, optional, default=True
        If True, normalize the molecule by checking aromaticity, adding
        explicit hydrogens, and renaming by IUPAC name.
    keep_confs : int, optional, default=None
        If None, apply the charges to the provided conformation and return
        this conformation, unless no conformation is present.
        Otherwise, return some or all of the generated
        conformations. If -1, all generated conformations are returned.
        Otherwise, keep_confs = N will return an OEMol with up to N
        generated conformations.  Multiple conformations are still used to
        *determine* the charges.
    Returns
    -------
    charged_copy : OEMol
        A molecule with OpenEye's recommended AM1BCC charge selection scheme.
    Notes
    -----
    Roughly follows
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    """

    # If there is no geometry, return at least one conformation.
    if molecule.GetConfs() == 0:
        keep_confs = 1

    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oequacpac = import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))

    if normalize:
        molecule = normalize_molecule(molecule)
    else:
        molecule = oechem.OEMol(molecule)

    charged_copy = generate_conformers(molecule, max_confs=max_confs, strictStereo=strictStereo)  # Generate up to max_confs conformers

    # 2017.2.1 Release new charging function
    status = oequacpac.OEAssignCharges(charged_copy, oequacpac.OEAM1BCCCharges())

    if not status:
        raise(RuntimeError("OEAssignCharges returned error code %d" % status))

    #Determine conformations to return
    if keep_confs == None:
        #If returning original conformation
        original = molecule.GetCoords()
        #Delete conformers over 1
        for k, conf in enumerate( charged_copy.GetConfs() ):
            if k > 0:
                charged_copy.DeleteConf(conf)
        #Copy coordinates to single conformer
        charged_copy.SetCoords( original )
    elif keep_confs > 0:
        #Otherwise if a number is provided, return this many confs if available
        for k, conf in enumerate( charged_copy.GetConfs() ):
            if k > keep_confs - 1:
                charged_copy.DeleteConf(conf)
    elif keep_confs == -1:
        #If we want all conformations, continue
        pass
    else:
        #Not a valid option to keep_confs
        raise(ValueError('Not a valid option to keep_confs in get_charges.'))

    return charged_copy

class GenMolStructure(object):
    """
    Generates parameterized ParmEd structure of the molecule with a chosen forcefield

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The openeye molecule to be parameterized
    forcefield : str
        String specifying the forcefield parameters to be used.

    Returns
    ---------
    packedmol : openeye.oechem.OEMol
        Openeye molecule with the ParmEd Structure attached.
    """

    def __init__(self, molecule, forcefield):
        if not forcefield in ['SMIRNOFF', 'GAFF', 'GAFF2']:
            raise RuntimeError('Selected forcefield %s is not GAFF/GAFF2/SMIRNOFF' % forcefield)
        else:
            self.molecule = molecule
            self.forcefield = str(forcefield).strip()
            self.structure = None

    def checkTleap(self):
        #Try to check if tleap is going to fail
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

    def smirnoffStructure(self, molecule=None):
        from smarty.forcefield import ForceField
        from smarty.forcefield_utils import create_system_from_molecule
        if not molecule:
            molecule = self.molecule

        try:
            ff = get_data_filename('smirff99Frosst','smirff99Frosst.ffxml')
            with open(ff) as ffxml:
                mol_ff = ForceField(ffxml)
        except:
            raise RuntimeError('Error opening {}'.format(ff))

        self.checkCharges(molecule)
        mol_top, mol_sys, mol_pos = create_system_from_molecule(mol_ff, molecule)
        molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)

        return molecule_structure

    def gaffStructure(self, molecule=None, forcefield=None):
        if not molecule:
            molecule = self.molecule
        if not forcefield:
            forcefield = self.forcefield

        #Try to check if tleap is going to fail
        self.checkTleap()
        self.checkCharges(molecule)

        # Determine formal charge (antechamber needs as argument)
        chg = 0
        for atom in molecule.GetAtoms():
            chg+=atom.GetFormalCharge()

        # Write out mol to a mol2 file to process via AmberTools
        mol2file = tempfile.NamedTemporaryFile(suffix='.mol2')
        mol2filename = mol2file.name
        with oechem.oemolostream(mol2filename) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, molecule)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing molecule %s to mol2." % molecule.GetTitle())

        # Run antechamber to type and parmchk for frcmod
        # requires openmoltools 0.7.5 or later, which should be conda-installable via omnia
        gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber('ligand', mol2filename,
                                                                                 gaff_version=forcefield.lower(),
                                                                                 net_charge=chg)

        # Run tleap using specified forcefield
        prmtop, inpcrd = openmoltools.amber.run_tleap('ligand', gaff_mol2_filename,
                                                      frcmod_filename,
                                                      leaprc = 'leaprc.%s' % forcefield.lower())

        # Load via ParmEd
        molecule_structure = parmed.amber.AmberParm( prmtop, inpcrd )

        return molecule_structure

    def parameterize(self):
        if self.forcefield == 'SMIRNOFF':
            structure = self.smirnoffStructure()
        elif self.forcefield in ['GAFF', 'GAFF2']:
            structure = self.gaffStructure()
        self.structure = structure
        return self.structure