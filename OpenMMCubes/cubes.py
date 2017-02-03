import time
import traceback
import numpy as np
from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES
from OpenMMCubes.ports import OpenMMSystemOutput, OpenMMSystemInput
from simtk import unit, openmm
from simtk.openmm import app

from openeye import oechem
import os, smarty, parmed, pdbfixer
from openmoltools import forcefield_generators
from simtk.openmm import XmlSerializer
from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
import OpenMMCubes.utils as utils

#For parallel, import and inherit from ParallelOEMolComputeCube
class OpenMMComplexSetup(OEMolComputeCube):
    title = "OpenMMComplexSetup"
    description = """
    Set up protein:ligand complex for simulation with OpenMM.

    This cube will generate an OpenMM System object containing
    a TIP3P solvated protein:ligand complex. The complex
    will be stored into a complex.oeb.gz file and streamed into the OpenMMSimulation cube.
    """
    classification = [
        ["OpenMM", "ProtLigComplex Setup"],
    ]
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
        default=5,
        help_text="Padding around protein for solvent box (angstroms)",
    )

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=50,
        help_text="Salt concentration (millimolar)",
    )

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
    #    required=True,
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein'
    )

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
    #   required=True,
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent'
    )

    def begin(self):
        pdbfilename = 'protein.pdb'
        # Write the protein to a PDB
        if in_orion():
            stream = StreamingDataset(self.args.protein, input_format=".pdb")
            stream.download_to_file(pdbfilename)
        else:
            protein = oechem.OEMol()
            with oechem.oemolistream(self.args.protein) as ifs:
                if not oechem.OEReadMolecule(ifs, protein):
                    raise RuntimeError("Error reading molecule")

        # Read the PDB file into an OpenMM PDBFile object
        self.proteinpdb = app.PDBFile(self.args.protein)

    def check_tagdata(self, mol):
        # Ensure tagged generic data is retained across cubes
        if 'idtag' not in mol.GetData().keys():
            raise RuntimeError('Could not find idtag for molecule')
        else:
            idtag =  mol.GetData(oechem.OEGetTag('idtag'))
            self.idtag = idtag
        # Regenerate Parmed Molecule structure
        if 'system' not in mol.GetData().keys():
            raise RuntimeError("Could not find system for molecule")
        else:
            intake = OpenMMSystemInput('intake')
            systag = oechem.OEGetTag('system')
            system = intake.decode(mol.GetData(systag))
            positions = utils.getPositionsFromOEMol(mol)
            topology = utils.generateTopologyFromOEMol(mol)
            structure = parmed.openmm.load_topology(topology,
                                                    system,
                                                    xyz=positions)
            structure.residues[0].name = "MOL"
            self.structure = structure

        if not any([self.idtag, self.structure]):
            raise RuntimeError('Missing tagged generic data')
        else:
            return True

    def process(self, mol, port):
        if self.check_tagdata(mol):
            idtag = self.idtag
            self.outfname = 'output/{}-complex'.format(idtag)
            outfname = self.outfname
            molecule_structure = self.structure
        try:
            #Generate protein Structure object
            forcefield = app.ForceField(self.args.protein_forcefield, self.args.solvent_forcefield)
            protein_system = forcefield.createSystem( self.proteinpdb.topology )
            protein_structure = parmed.openmm.load_topology( self.proteinpdb.topology,
                                                             protein_system,
                                                             xyz=self.proteinpdb.positions )

            # Merge structures to prevent adding solvent in pocket
            pl_structure = protein_structure + molecule_structure
            self.log.info('{}-complex: {}'.format(idtag, pl_structure))

            # Retain positions and save
            pl_structure.positions = utils.combinePostions(protein_structure.positions,
                                            molecule_structure.positions)
            pl_structure.save(outfname+'-pl.tmp',format='pdb',overwrite=True)

            # Solvate with PDBFixer
            self.log.info('PDBFixer solvating {}-complex:'.format(idtag))
            self.log.info('\tpH = {}'.format(self.args.pH))
            self.log.info('\tpadding = {}'.format(unit.Quantity(self.args.solvent_padding, unit.angstroms)))
            self.log.info('\tionicStrength = {}'.format(unit.Quantity(self.args.salt_concentration, unit.millimolar)))
            fixer = pdbfixer.PDBFixer(outfname+'-pl.tmp')
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            if fixer.missingAtoms:
                self.log.info('Adding missing atoms: {}'.format(fixer.missingAtoms))
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
            tmp.strip(":MOL")
            tmp.save(outfname+'-nomol.tmp',format='pdb',overwrite=True)
            # Reload PDBFile
            nomol = app.PDBFile(outfname+'-nomol.tmp')
            # Regenerate openMM System to parameterize solvent
            nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)
            # Regenerate parameterized solvated protein structure
            solv_structure = parmed.openmm.load_topology(nomol.topology,
                                                        nomol_system,
                                                        xyz=nomol.positions,
                                                        box=full_box)

            # Remerge with ligand structure
            full_structure = solv_structure + molecule_structure
            # Restore box dimensions
            full_structure.box = full_box
            # Save full structure
            full_structure.save(outfname+'.pdb', overwrite=True)
            self.log.info('Solvated {}-complex {}'.format(idtag, full_structure))
            self.log.info('\tBox = {}'.format(full_structure.box))

            # Regenerate OpenMM system with parmed
            system = full_structure.createSystem(nonbondedMethod=app.PME,
                                                nonbondedCutoff=10.0*unit.angstroms,
                                                constraints=app.HBonds)

            # Pack solvated complex into oeb and emit system
            complex_mol = oechem.OEMol()
            output = OpenMMSystemOutput('output')
            with oechem.oemolistream(outfname+'.pdb') as ifs:
                if not oechem.OEReadMolecule(ifs, complex_mol):
                    raise RuntimeError("Error reading {}.pdb".format(outfname))

                with oechem.oemolostream(outfname+'.oeb.gz') as ofs:
                    complex_mol.SetData(oechem.OEGetTag('idtag'), idtag)
                    complex_mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                    res = oechem.OEWriteConstMolecule(ofs, complex_mol)
                    if res != oechem.OEWriteMolReturnCode_Success:
                        raise RuntimeError("Error writing {}.oeb.gz".format(outfname))
                    else:
                        self.log.info('Saved System to: {}.oeb.gz'.format(outfname))
            self.success.emit(complex_mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

    def end(self):
        #Clean up
        os.remove(self.outfname+'-pl.tmp')
        os.remove(self.outfname+'-nomol.tmp')

class OpenMMSimulation(OEMolComputeCube):
    title = "Run simulation in OpenMM"
    description = """
    Run simulation with OpenMM for protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the protein:ligand complex, reconstruct the OpenMM System object,
    minimize the system and run 1000 MD steps at 300K.
    The potential energies are evaluated every 100 steps and stored to a log file.
    Stdout is a progress/benchmark timings reporter.
    The OpenMM System, State, and log file are attached to the OEMol and saved
    to the file simulation.oeb.gz.

    The simulation.oeb.gz file, containing the State can then be reused to
    restart the MD simulation.
    """
    classification = [
        ["Testing", "OpenMM"],
        ["Testing", "Simulation"],
    ]
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
        default=1000,
        help_text="Number of MD steps")

    complex_mol = parameter.DataSetInputParameter(
        'complex_mol',
        #default='complex.oeb.gz',
        help_text='Single protein to Dock Against')

    def __init__(self, *args, **kwargs):
        super(OpenMMSimulation, self).__init__(*args, **kwargs)
        self._setup = False

    def begin(self):
        if self._setup:
            return
        self._setup = True

    def check_tagdata(self, mol):
        # Ensure tagged generic data is retained across cubes
        if 'idtag' not in mol.GetData().keys():
            raise RuntimeError('Could not find idtag for molecule')
        else:
            idtag =  mol.GetData(oechem.OEGetTag('idtag'))
            self.idtag = idtag
        # Regenerate Parmed Molecule structure
        if 'system' not in mol.GetData().keys():
            raise RuntimeError("Could not find system for molecule")
        else:
            intake = OpenMMSystemInput('intake')
            systag = oechem.OEGetTag('system')
            system = intake.decode(mol.GetData(systag))
            positions = utils.getPositionsFromOEMol(mol)
            topology = utils.generateTopologyFromOEMol(mol)
            self.system = system
            self.positions = positions
            self.topology = topology

        # Check if mol has State data attached
        if 'state' in mol.GetData().keys():
            self.log.info('Found a saved State, restarting simulation')
            mol.GetData(oechem.OEGetTag('state'))
            serialized_state = mol.GetData(oechem.OEGetTag('state'))
            state = openmm.XmlSerializer.deserialize( serialized_state )
            self.state = state
        else:
            self.state = None

        if not any([self.idtag, self.system, self.positions, self.topology]):
            raise RuntimeError('Missing tagged generic data')
        else:
            return True

    def setReporters(self):
        from sys import stdout
        progress_reporter = app.StateDataReporter(stdout,
                                            reportInterval=100, totalSteps=1000,
                                            time=True, speed=True, progress=True,
                                            elapsedTime=True, remainingTime=True)

        state_reporter = app.StateDataReporter(self.outfname+'.log',
                                            reportInterval=100, step=True,
                                            potentialEnergy=True, totalEnergy=True,
                                            volume=True, temperature=True)
        self.reporters = [progress_reporter, state_reporter]
        return self.reporters

    def process(self, complex_mol, port):
        if self.check_tagdata(complex_mol):
            idtag = self.idtag
            self.outfname = 'output/{}-simulation'.format(idtag)
            outfname = self.outfname
            system = self.system
            positions = self.positions
            topology = self.topology
        try:
            # Initialize Simulation
            integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
            simulation = app.Simulation(topology, system, integrator)
            #simulation = app.Simulation(topology, system, self.integrator, openmm.Platform.getPlatformByName('CPU'))
            platform = simulation.context.getPlatform().getName()
            self.log.info('Running OpenMMSimulation on Platform {}'.format(platform))

            # Check if mol has State data attached
            if self.state:
                simulation.context.setState(self.state)
                self.outfname = 'output/{}-restart'.format(idtag)
                outfname = self.outfname
            else:
                # Set initial positions and velocities then minimize
                simulation.context.setPositions(positions)
                simulation.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)
                init = simulation.context.getState(getEnergy=True)
                self.log.info('Initial energy is {}'.format(init.getPotentialEnergy()))
                # Temporarily, place some restrictions on minization to run faster
                self.log.info('Minimizing {} system...'.format(idtag))
                if in_orion():
                    simulation.minimizeEnergy(tolerance=unit.Quantity(10.0,unit.kilojoules/unit.moles),maxIterations=100)
                else:
                    simulation.minimizeEnergy()
                st = simulation.context.getState(getPositions=True,getEnergy=True)
                self.log.info('Minimized energy is {}'.format(st.getPotentialEnergy()))


            #Append Reporters to simulation
            reporters = self.setReporters()
            for rep in reporters:
                simulation.reporters.append(rep)

            self.log.info('Running {} MD steps at {}K'.format(self.args.steps, self.args.temperature))

            simulation.step(self.args.steps)
            outlog = open(outfname+'.log', 'r')
            self.log.info(outlog.read())

            # Save serialized State object
            state = simulation.context.getState( getPositions=True,
                                              getVelocities=True,
                                              getParameters=True )

            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            self.log.info('Saving to {}'.format(outfname+'.oeb.gz'))
            with oechem.oemolostream(outfname+'.oeb.gz') as ofs:
                complex_mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                complex_mol.AddData(oechem.OEGetTag('state'), output.encode(state))
                complex_mol.AddData(oechem.OEGetTag('log'), outlog.read())
                oechem.OEWriteConstMolecule(ofs, complex_mol)
            self.success.emit(complex_mol)
            outlog.close()

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                complex_mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(complex_mol)
