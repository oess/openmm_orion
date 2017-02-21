import io, os, time, traceback, base64, smarty, parmed, pdbfixer, mdtraj, pickle
from openeye import oechem
from sys import stdout
import numpy as np
from simtk import unit, openmm
from simtk.openmm import app

from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort, ParallelOEMolComputeCube
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES

from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
import OpenMMCubes.utils as utils
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )
from OpenMMCubes.utils import download_dataset_to_file, get_data_filename


class OpenMMComplexSetup(OEMolComputeCube):
    title = "OpenMMComplexSetup"
    description = """
    Set up protein:ligand complex for simulation with OpenMM.

    This cube will generate an OpenMM System containing
    a TIP3P solvated protein:ligand complex. The complex
    will be stored into a <idtag>-complex.oeb.gz file, with the System and Structure
    attached and streamed into the OpenMMSimulation cube.
    """
    classification = ["Complex Setup"]
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    protein = parameter.DataSetInputParameter(
        'protein',
        required=True,
        help_text='Protein PDB file')

    pH = parameter.DecimalParameter(
        'pH',
        default=7.0,
        help_text="Solvent pH used to select appropriate protein protonation state.",
    )

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10,
        help_text="Padding around protein for solvent box (angstroms)",
    )

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=50,
        help_text="Salt concentration (millimolar)",
    )

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein'
    )

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent'
    )

    def begin(self):
        pdbfilename = 'protein.pdb'
        protein = oechem.OEMol()
        self.args.protein = download_dataset_to_file(self.args.protein)
        with oechem.oemolistream(self.args.protein) as ifs:
            if not oechem.OEReadMolecule(ifs, protein):
                raise RuntimeError("Error reading protein")
        with oechem.oemolostream(pdbfilename) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, protein)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing protein: {}".format(res))

        # Read the PDB file into an OpenMM PDBFile object
        self.proteinpdb = app.PDBFile(pdbfilename)

    def process(self, mol, port):
        try:
            req_tags = ['idtag', 'structure']
            if utils.OEPackMol.checkTags(mol, req_tags):
                gd = utils.OEPackMol.unpack(mol)
                outfname = '{}-complex'.format(gd['idtag'])

            #Generate protein Structure object
            forcefield = app.ForceField(self.args.protein_forcefield, self.args.solvent_forcefield)
            protein_system = forcefield.createSystem( self.proteinpdb.topology )
            protein_structure = parmed.openmm.load_topology( self.proteinpdb.topology,
                                                             protein_system,
                                                             xyz=self.proteinpdb.positions )

            # Merge structures to prevent adding solvent in pocket
            pl_structure = protein_structure + gd['structure']
            self.log.info('{}-complex: {}'.format(gd['idtag'], pl_structure))

            # Retain positions and save
            pl_structure.positions = utils.combinePostions(protein_structure.positions,
                                            gd['structure'].positions)
            pl_structure.save(outfname+'-pl.tmp',format='pdb',overwrite=True)

            # Solvate with PDBFixer
            self.log.info('PDBFixer solvating {}-complex:'.format(gd['idtag']))
            self.log.info('\tpH = {}'.format(self.args.pH))
            self.log.info('\tpadding = {}'.format(unit.Quantity(self.args.solvent_padding, unit.angstroms)))
            self.log.info('\tionicStrength = {}'.format(unit.Quantity(self.args.salt_concentration, unit.millimolar)))
            fixer = pdbfixer.PDBFixer(outfname+'-pl.tmp')
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.replaceNonstandardResidues()
            #fixer.removeHeterogens(False)
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(self.args.pH)
            fixer.addSolvent(padding=unit.Quantity(self.args.solvent_padding, unit.angstroms),
                            ionicStrength=unit.Quantity(self.args.salt_concentration, unit.millimolar)
                            )

            # Load PDBFixer object back to Structure
            tmp = parmed.openmm.load_topology(fixer.topology, xyz=fixer.positions)
            #Store positions, topology, and box vectors for solvated system
            full_positions = tmp.positions
            full_topology = tmp.topology
            full_box = tmp.box
            # Remove ligand from protein Structure by AmberMask selection
            tmp.strip(":LIG")
            tmp.save(outfname+'-nomol.tmp',format='pdb',overwrite=True)
            # Reload PDBFile
            nomol = app.PDBFile(outfname+'-nomol.tmp')
            nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)
            # Regenerate parameterized solvated protein structure
            solv_structure = parmed.openmm.load_topology(nomol.topology,
                                                        nomol_system,
                                                        xyz=nomol.positions,
                                                        box=full_box)

            # Remerge with ligand structure
            full_structure = solv_structure + gd['structure']
            # Restore box dimensions
            full_structure.box = full_box
            self.log.info('Solvated {}-complex {}'.format(gd['idtag'], full_structure))
            self.log.info('\tBox = {}'.format(full_structure.box))

            packedmol = utils.OEPackMol.pack(mol, full_structure)
            self.success.emit(packedmol)

            #Cleanup
            utils.cleanup(['protein.pdb', outfname+'-pl.tmp',
                       outfname+'-nomol.tmp'])

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

class OpenMMSimulation(OEMolComputeCube):
    title = "Run simulation in OpenMM"
    description = """
    Run simulation with OpenMM for protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the protein:ligand complex, reconstruct the OpenMM System,
    minimize the system, save the minimized PDB, and run 1000 MD steps at 300K.
    The potential energies are evaluated every 1000 steps and stored to a log file.
    Stdout is a progress/benchmark timings reporter every 1000 steps.
    The Structure, OpenMM System, State, and log file are attached to the OEMol and
    saved to the file simulation.oeb.gz.

    The simulation.oeb.gz file, containing the State can then be reused to
    restart the MD simulation.
    """
    classification = ["Simulation"]
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)"
    )
    steps = parameter.IntegerParameter(
        'steps',
        default=50000,
        help_text="Number of MD steps")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=1000,
        help_text="Step interval for reporting data."
    )

    def begin(self):
        pass

    def process(self, mol, port):
        try:
            req_tags = ['idtag', 'structure']
            if utils.OEPackMol.checkTags(mol, req_tags):
                gd = utils.OEPackMol.unpack(mol)
                outfname = '{}-simulation'.format(gd['idtag'])

            simulation = utils.genSimFromStruct(gd['structure'], self.args.temperature)
            # Check if mol has State data attached
            if 'state' in gd.keys():
                simulation.context.setState(gd['state'])
            else:
                simulation = utils.minimizeSimulation(simulation)

            reporters = utils.setReporters(self.args.reporter_interval, self.args.steps,
                                           outfname)
            #Append Reporters to simulation
            for rep in reporters:
                simulation.reporters.append(rep)

            self.log.info('Running {} MD steps at {}K'.format(self.args.steps, self.args.temperature))
            simulation.step(self.args.steps)

            packedmol = utils.OEPackMol.pack(mol, simulation)
            self.success.emit(packedmol)
            tmpfiles = [ gd['idtag']+'-simulation.log', gd['idtag']+'-simulation.nc' ]
            utils.cleanup(tmpfiles)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(mol)
