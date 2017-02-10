import io, os, traceback, string, random, parmed
from openeye import oechem, oedocking, oeomega
from floe.api import (
    parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube, MoleculeInputPort,
    StringParameter, MoleculeOutputPort
)
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES

from LigPrepCubes.ports import (
    CustomMoleculeInputPort, CustomMoleculeOutputPort)
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )

from smarty.forcefield import ForceField
from smarty.forcefield_utils import create_system_from_molecule
from OpenMMCubes.utils import download_dataset_to_file, get_data_filename

def _generateRandomID(size=5, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class SetIDTagfromTitle(OEMolComputeCube):
    title = "SetIDTagfromTitle"
    description = """
    Attach IDname to OEMol tag.
    """
    classification = [["Testing", "Ligand Preparation"]]
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
    classification = [["Testing", "Ligand Preparation"]]
    tags = [tag for lists in classification for tag in lists]

    molecule_forcefield = parameter.DataSetInputParameter(
        'molecule_forcefield',
        default='smirff99Frosst.ffxml',
        help_text='Forcefield FFXML file for molecule')

    def begin(self):
        try:
            dfname = get_data_filename(self.args.molecule_forcefield)
            with open(dfname) as ffxml:
                self.mol_ff = ForceField(ffxml)
        except:
            raise RuntimeError('Error opening {}'.format(dfname))

    def process(self, mol, port):
        # Create a copy incase of error
        init_mol = oechem.OEMol(mol)
        try:
            mol_top, mol_sys, mol_pos = create_system_from_molecule(self.mol_ff, mol)
            molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)
            molecule_structure.residues[0].name = "MOL"

            # Encode System/Structure, Attach to mol
            sys_out = OpenMMSystemOutput('sys_put')
            struct_out = ParmEdStructureOutput('struct_out')
            mol.SetData(oechem.OEGetTag('system'), sys_out.encode(mol_sys))
            mol.SetData(oechem.OEGetTag('structure'), struct_out.encode(molecule_structure))
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
    #classification = [["Testing", "Output"]]
    title = "OEBSinkCube"
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
        if 'idtag' in mol.GetData().keys():
            idtag = mol.GetData(oechem.OEGetTag('idtag'))
        outfname = '{}/{}-{}.oeb.gz'.format(self.args.directory,
                                           idtag, self.args.suffix)
        with oechem.oemolostream(outfname) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, mol)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing {}.oeb.gz".format(outfname))
            else:
                self.log.info('Saving to {}'.format(outfname))
