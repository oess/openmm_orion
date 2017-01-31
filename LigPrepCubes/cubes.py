import traceback
from floe.api import (
    parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube, MoleculeInputPort,
    StringParameter, MoleculeOutputPort
)
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES
from openeye import oechem, oedocking, oeomega
from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort

class Attachffxml(ParallelOEMolComputeCube):
    title = "Attach FFXML to OE molecules"
    description = """
    Attach FFXML to OE molecules
    """
    classification = [
        ["OpenEye", "Ligand Preparation"],
    ]
    tags = [tag for lists in classification for tag in lists]

    molecule_forcefield = parameter.DataSetInputParameter(
        'molecule_forcefield',
        required=True,
        help_text='Forcefield FFXML file for molecule')

    def begin(self):        # Write the protein to a PDB
        if in_orion():
            pass
        else:
            self.ffxml = open(self.args.molecule_forcefield, 'rb')

    def process(self, mol, port):
        cubename = '[{}]'.format( str(self.name) )
        try:
            mol.SetData(oechem.OEGetTag('forcefield'), self.ffxml.read())
            self.success.emit(mol)
        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

class FREDDocking(ParallelOEMolComputeCube):
    title = "Docking Molecules"
    description = """
    Dock OE molecules
    """
    classification = [
        ["OpenEye", "Ligand Preparation"],
    ]
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    receptor = parameter.DataSetInputParameter(
        'receptor',
        required=True,
        help_text='Receptor OEB File')

    def __init__(self, *args, **kwargs):
        super(FREDDocking, self).__init__(*args, **kwargs)
        self._setup = False

    def begin(self):
        if self._setup:
            return
        #Read in the Receptor
        receptor = oechem.OEGraphMol()
        #apofname = 'epoxide_hydrolase_apo_receptor.oeb.gz'
        if not oedocking.OEReadReceptorFile(receptor, self.args.receptor):
            OEThrow.Fatal("Unable to read receptor")

        #Initialize Docking
        dock_method = oedocking.OEDockMethod_Hybrid
        if not oedocking.OEReceptorHasBoundLigand(receptor):
            oechem.OEThrow.Warning("No bound ligand, switching OEDockMethod to ChemGauss4.")
            dock_method = oedocking.OEDockMethod_Chemgauss4
        dock_resolution = oedocking.OESearchResolution_Default
        self.sdtag = oedocking.OEDockMethodGetName(dock_method)
        self.dock = oedocking.OEDock(dock_method, dock_resolution)
        if not self.dock.Initialize(receptor):
            raise Exception("Unable to initialize Docking with {0}".format(self.args.receptor))
        self._setup = True

    def clean(self, mol):
        mol.DeleteData('CLASH')
        mol.DeleteData('CLASHTYPE')
        mol.GetActive().DeleteData('CLASH')
        mol.GetActive().DeleteData('CLASHTYPE')

    def process(self, mol, port):
        molname = mol.GetTitle()

        #Set AtomTypes
        oechem.OETriposAtomNames(mol)
        oechem.OETriposAtomTypeNames(mol)

        cubename = '[{}]'.format( str(self.name) )
        dockedMol = oechem.OEMol()
        self.log.info("{} has {} conformers".format(molname, mol.NumConfs()))
        res = self.dock.DockMultiConformerMolecule(dockedMol, mol)
        if res == oedocking.OEDockingReturnCode_Success:

            oedocking.OESetSDScore(dockedMol, self.dock, self.sdtag)
            if oechem.OEHasSDData(dockedMol, self.sdtag):
                method = str(oechem.OEGetSDData(dockedMol, self.sdtag))
                oechem.OEDeleteSDData(dockedMol, self.sdtag)
                oechem.OESetSDData(dockedMol, self.sdtag, "%s".format(method))
            oechem.OESetSDData(dockedMol, "__ParentMol", self.sdtag)

            self.dock.AnnotatePose(dockedMol)
            self.log.info("%s score = %f" % (molname, self.dock.ScoreLigand(dockedMol)))

            self.clean(dockedMol)
            self.success.emit(dockedMol)

        else:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)
        return dockedMol
