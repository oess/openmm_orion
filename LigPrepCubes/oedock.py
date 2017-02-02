import traceback
from floe.api import (
    parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube, MoleculeInputPort,
    StringParameter, MoleculeOutputPort
)
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES
from openeye import oechem, oedocking, oeomega
from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort

class FREDDocking(OEMolComputeCube):
    title = "FREDDocking"
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

    def begin(self):
        #Read in the Receptor
        receptor = oechem.OEGraphMol()
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

    def clean(self, mol):
        mol.DeleteData('CLASH')
        mol.DeleteData('CLASHTYPE')
        mol.GetActive().DeleteData('CLASH')
        mol.GetActive().DeleteData('CLASHTYPE')

    def process(self, mcmol, port):
        dockedMol = oechem.OEMol()
        res = self.dock.DockMultiConformerMolecule(dockedMol, mcmol)
        if res == oedocking.OEDockingReturnCode_Success:

            oedocking.OESetSDScore(dockedMol, self.dock, self.sdtag)
            if oechem.OEHasSDData(dockedMol, self.sdtag):
                method = str(oechem.OEGetSDData(dockedMol, self.sdtag))
                oechem.OEDeleteSDData(dockedMol, self.sdtag)
                oechem.OESetSDData(dockedMol, self.sdtag, "%s".format(method))
            oechem.OESetSDData(dockedMol, "__ParentMol", self.sdtag)

            # Get top scoring pose
            self.dock.AnnotatePose(dockedMol)
            self.log.info("FRED %s score = %f" % (dockedMol.GetTitle(), self.dock.ScoreLigand(dockedMol)))
            self.clean(dockedMol)
            self.success.emit(dockedMol)

        else:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)
        return dockedMol

    def end(self):
        pass
