import traceback, string, random, os, io
from floe.api import (
    parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube, MoleculeInputPort,
    StringParameter, MoleculeOutputPort
)
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES
from openeye import oechem, oedocking, oeomega
from LigPrepCubes.ports import (
    CustomMoleculeInputPort, CustomMoleculeOutputPort,
    ParmEdStructureInput, ParmEdStructureOutput)
from OpenMMCubes.ports import OpenMMSystemOutput, OpenMMSystemInput

def _generateRandomID(size=5, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class SetIDTagfromTitle(OEMolComputeCube):
    title = "SetIDTagfromTitle"
    description = """
    Attach IDname to OEMol tag.
    """
    classification = [["OpenEye", "Ligand Preparation"]]
    tags = [tag for lists in classification for tag in lists]

    def process(self, mol, port):
        #Check for OEMol title for ID labeling
        if not mol.GetTitle():
            idtag = _generateRandomID()
            oechem.OEThrow.Warning('No title found, setting to {}'.format(idtag))
            mol.SetTitle(idtag)
        else:
            idtag = mol.GetTitle()
        try:
            # Set AtomTypes
            oechem.OETriposAtomNames(mol)
            oechem.OETriposAtomTypeNames(mol)
            mol.AddData(oechem.OEGetTag('idtag'), mol.GetTitle())
            self.success.emit(mol)
        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)


class SMIRFFParameterization(OEMolComputeCube):
    title = "Attach FFXML to OE molecules"
    description = """
    Parameterize the ligand with the smirff99Frosst.ffxml parameters,
    which is parsed with smarty. Attach the System to the OEMol.
    """
    classification = [["OpenEye", "Ligand Preparation"]]
    tags = [tag for lists in classification for tag in lists]

    molecule_forcefield = parameter.DataSetInputParameter(
        'molecule_forcefield',
        required=True,
        help_text='Forcefield FFXML file for molecule')

    def begin(self):
        try:
            ffxml = open(self.args.molecule_forcefield, 'rb')
        except:
            raise RuntimeError('Error opening {}'.format(self.args.molecule_forcefield))
        ffxml.close()

    def process(self, mol, port):
        # Create a copy incase of error
        init_mol = oechem.OEMol(mol)
        try:
            from smarty.forcefield import ForceField
            from smarty.forcefield_utils import create_system_from_molecule
            with open( self.args.molecule_forcefield, 'r') as ffxml:
                mol_ff = ForceField( ffxml )
            mol_topology, mol_system, mol_positions = create_system_from_molecule(mol_ff, mol)

            output = OpenMMSystemOutput('output')
            mol.SetData(oechem.OEGetTag('system'), output.encode(mol_system))
            self.success.emit(mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            init_mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(init_mol)


class OEBSinkCube(SinkCube):
    """
    A custom sink cube that writes molecules to a oeb.gz
    """
    classification = [["Output"]]
    title = "Dataset Writer"
    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')

    directory = parameter.StringParameter('directory',
                                         default='output',
                                         description='Directory name')
    suffix = parameter.StringParameter('suffix',
                                        required=True,
                                        description='suffix to append')
    def begin(self):
        if not os.path.exists(self.args.directory):
            os.makedirs(self.args.directory)

    def write(self, mol, port):
        outfname = '{}/{}-{}.oeb.gz'.format(self.args.directory,
                                           mol.GetTitle(), self.args.suffix)
        with oechem.oemolostream(outfname) as ofs:
            oechem.OEWriteConstMolecule(ofs, mol)
        self.log.info('Saving to {}'.format(outfname))
