import io, os, random, string, subprocess, tempfile, traceback
import openmoltools, parmed
from openeye import oechem, oedocking, oeomega
import OpenMMCubes.utils as utils
from LigPrepCubes import ff_utils
from floe.api import OEMolComputeCube, parameter

def _generateRandomID(size=5, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class ChargeMCMol(OEMolComputeCube):
    title = "Charge Multiconf. Molecules"
    version = "0.0.2"
    classification = [ ["Ligand Preparation", "OEChem", "Add Hydrogen"],
    ["Ligand Preparation", "OEChem", "Check Aromaticity"],
    ["Ligand Preparation", "OEChem", "IUPAC"],
    ["Ligand Preparation", "OMEGA", "Conformer Generation"],
    ["Ligand Preparation", "QUACPAC", "Charge Assignment"]]
    tags = ['Openmoltools', 'OMEGA', 'QUACPAC']
    description = """
    Calls openmoltools to perform the following:
    (1) 'normalize_molecule': checks aromaticity, add explicit hydrogens and renaming by IUPAC.
    (2) Generate multiple conformers with OMEGA.
    (3) Assigns partial charges with oequacpac.OEAssignCharges (req: OpenEye-toolkits: 2017.2.1).

    Input:
    -------
    oechem.OEMol - Streamed-in uncharged molecule with no hydrogens.

    Output:
    -------
    oechem.OEMCMol - Emits a charged multi-conformer molecule with attachments:
        - SDData Tags: { IUPAC : str, IDTag : str }
        - Generic Tags: { IDTag: str }
    """


    max_conf = parameter.IntegerParameter(
        'max_conf',
        default=800,
        help_text="Max number of conformers")


    keep_conf = parameter.IntegerParameter(
        'keep_conformers',
        default=None,
        help_text="Select the number of conformers to keep")

    legacy = parameter.BooleanParameter(
        'legacy',
        default=True,
        help_text="(Default=True) Use legacy charging engine OEAssignPartialCharges. False to use OEAssignCharges (Requires: OEToolkits 2017.2.1).")

    def begin(self):
        self.opt = vars(self.args)

    def process(self, mol, port):
        try:
            if not mol.GetTitle():
                idtag = _generateRandomID()
                self.log.warn('Mol title not found, setting to {}'.format(idtag))
            else:
                # Store the IDTag from the SMILES file.
                idtag = mol.GetTitle()
            import openmoltools
            charged_mol = openmoltools.openeye.get_charges(mol,
                            max_confs=self.opt['max_conf'],
                            strictStereo=True, normalize=True,
                            keep_confs=self.opt['keep_conf'],
                            legacy=False)
            #Generate the charged molecule, keeping the first conf.
            #charged_mol = ff_utils.assignCharges(mol,
            #                max_confs=self.opt['max_conformers'],
            #                strictStereo=True, normalize=True,
            #                keep_confs=self.opt['keep_conformers'])

            # Store the IUPAC name from normalize_molecule
            iupac = [ charged_mol.GetTitle().strip() ]
            # Pack as list incase of commas in IUPUC

            # Keep it as an SD Tag
            oechem.OESetSDData(charged_mol, 'IUPAC', str(iupac))
            # Reset the charged mol title to the original IDTag
            charged_mol.SetTitle(idtag)
            charged_mol.SetData(oechem.OEGetTag('IDTag'), idtag)
            oechem.OESetSDData(charged_mol, 'IDTag', idtag)

            self.success.emit(charged_mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

class LigandParameterization(OEMolComputeCube):
    title = "Ligand Parameterization"
    version = "0.0.3"
    classification = [ ["Ligand Preparation", "SMARTY", "Forcefield Assignment"],
    ["Ligand Preparation", "AMBER", "Forcefield Assignment"]]
    tags = ['Openmoltools', 'ParmEd', 'SMARTY', 'SMIRNOFF', 'GAFF']
    description = """
    Parameterize the ligand with the chosen forcefield.
    Supports GAFF/GAFF2/SMIRNOFF.
    Generate a parameterized parmed Structure of the molecule.

    Input:
    -------
    oechem.OEMol - Streamed-in charged molecule with explicit hydrogens.

    Output:
    -------
    oechem.OEMol - Emits molecule with attachments:
        - SDData Tags: { NumAtoms : str, FF : str, Structure: str <parmed.Structure> }
        - Generic Tags: { Structure : parmed.Structure (base64-encoded) }
    """

    molecule_forcefield = parameter.StringParameter(
        'molecule_forcefield',
        required=True,
        default='GAFF2',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Forcefield to parameterize the molecule')

    def begin(self):
        if self.args.molecule_forcefield in ['GAFF', 'GAFF2']:
            ff_utils.checkTleap(self.args.molecule_forcefield)

    def process(self, mol, port):
        try:
            pmd = ff_utils.ParamLigStructure(mol, self.args.molecule_forcefield)
            molecule_structure = pmd.parameterize()
            molecule_structure.residues[0].name = "LIG"
            self.log.info(str(molecule_structure))

            oechem.OESetSDData(mol, 'NumAtoms', str(mol.NumAtoms()))
            oechem.OESetSDData(mol, 'Structure', str(molecule_structure))
            oechem.OESetSDData(mol, 'FF', str(self.args.molecule_forcefield) )
            mol.SetData(oechem.OEGetTag('IDTag'), mol.GetTitle())
            packedmol = utils.PackageOEMol.pack(mol, molecule_structure)

            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)


class FREDDocking(OEMolComputeCube):
    title = "FRED Docking"
    version = "0.0.2"
    classification = [ ["Ligand Preparation", "OEDock", "FRED"],
    ["Ligand Preparation", "OEDock", "ChemGauss4"]]
    tags = ['OEDock', 'FRED']
    description = """
    Dock molecules using the FRED docking engine against a prepared receptor file.
    Return the top scoring pose.

    Input:
    -------
    receptor - Requires a prepared receptor (oeb.gz) file of the protein to dock molecules against.
    oechem.OEMCMol - Expects a charged multi-conformer molecule on input port.

    Output:
    -------
    oechem.OEMol - Emits the top scoring pose of the molecule with attachments:
        - SDData Tags: { ChemGauss4 : pose score }
    """

    receptor = parameter.DataSetInputParameter(
        'receptor',
        required=True,
        help_text='Receptor OEB File')

    def begin(self):
        receptor = oechem.OEGraphMol()
        self.args.receptor = utils.download_dataset_to_file(self.args.receptor)
        if not oedocking.OEReadReceptorFile(receptor, str(self.args.receptor)):
            raise Exception("Unable to read receptor from {0}".format(self.args.receptor))

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
        try:
            dockedMol = oechem.OEMol()
            res = self.dock.DockMultiConformerMolecule(dockedMol, mcmol)

            if res == oedocking.OEDockingReturnCode_Success:
                oedocking.OESetSDScore(dockedMol, self.dock, self.sdtag)
                self.dock.AnnotatePose(dockedMol)
                score = self.dock.ScoreLigand(dockedMol)
                self.log.info("{} {} score = {:.4f}".format(self.sdtag, dockedMol.GetTitle(), score))
                oechem.OESetSDData(dockedMol, self.sdtag, "{}".format(score))
                self.clean(dockedMol)
                self.success.emit(dockedMol)
            else:
                self.log.error('{} - {}'.format(res, oedocking.OEDockingReturnCodeGetName(res)) )

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mcmol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mcmol)

    def end(self):
        pass
