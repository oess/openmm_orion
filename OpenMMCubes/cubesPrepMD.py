import io, os, time, traceback
from openeye import oechem
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

import json
from OpenMMCubes import plmd

class OpenMMmakePLmaskCube(OEMolComputeCube):
    title = "Generate the Protein-Ligand Mask used for MD restraints"
    description = """
    Generate the restraint mask for use with OpenMM for a protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and generate a dictionary to be used
    by OpenMM in placing restraints on various components of the system.
    The protein-ligand mask dictionary is attached to the OEMol and
    saved to the file PLmask.oeb.gz.

    Input parameter:
    ActSiteResNumSDTag (string): the SD Tag for the whitespace-delimited protein
    active site residue numbers (integers). These residue numbers must correspond
    to residue numbers in the streamed complex.oeb.gz file.
    """
    classification = ['MDPrep']
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    ActSiteResNumSDTag = parameter.StringParameter(
        'ActSiteResNumSDTag',
        default='ActiveSiteResNums',
        help_text="whitespace delimited list of integers corresponding to residue numbers")

    def process(self, complex_mol, port):
        try:
            # generate the protein-ligand mask (a python dict)
            atomPLmask = plmd.protLigMask(complex_mol, self.args.ActSiteResNumSDTag)
            # package the mask as a json object then attach to the molecule; emit
            dataTagForOpenMM_PLmask = 'OpenMM_PLmaskDict_json'
            jsonPLmask = json.dumps(atomPLmask, ensure_ascii=True)
            # attach the json object to the molecule; emit
            complex_mol.SetStringData( dataTagForOpenMM_PLmask, jsonPLmask)
            self.success.emit(complex_mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            complex_mol.SetData('OpenMMmakePLmask_Error', str(e))
            # Return failed mol
            self.failure.emit(complex_mol)

class OpenMMminimizeCube(OEMolComputeCube):
    title = 'Minimize with OpenMM using strong restraints'
    description = """
    Minimize the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and minimize it using
    strong restraints on the protein and ligand.

    Input parameters:
    steps (integer): the number of steps of minimization to apply.

    Input parameters:
    PLMaskSDTag (string): the SD Tag for the protein-ligand mask used
    to apply restraints to the minimization.
    """
    classification = ['PrepMDminimize']
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    steps = parameter.IntegerParameter(
        'steps',
        default=100,
        help_text="Number of MD steps")

    PLMaskSDTag = parameter.StringParameter(
        'PLMaskSDTag',
        default='OpenMM_PLmaskDict_json',
        help_text="SD Tag for the Protein-Ligand Mask for restraints")

    def check_tagdata(self, mol):
        if 'idtag' not in mol.GetData().keys():
            raise RuntimeError('Could not find idtag for molecule')
        else:
            idtag =  mol.GetData(oechem.OEGetTag('idtag'))
            self.idtag = idtag
        if 'system' not in mol.GetData().keys():
            raise RuntimeError("Could not find system for molecule")
        else:
            sys_in = OpenMMSystemInput('sys_in')
            sys_tag = oechem.OEGetTag('system')
            system = sys_in.decode(mol.GetData(sys_tag))
            self.system = system
        if 'structure' not in mol.GetData().keys():
            raise RuntimeError('Could not find structure for molecule')
        else:
            struct_in = ParmEdStructureInput('struct_in')
            struct_tag = oechem.OEGetTag('structure')
            structure = struct_in.decode(mol.GetData(struct_tag))
            self.structure = structure

        if not any([self.idtag, self.system, self.structure]):
            raise RuntimeError('Missing tagged generic data')
        else:
            return True

    def begin(self):
        if not os.path.exists('./output'):
            os.makedirs('./output')
        return

    def process(self, complex_mol, port):
        try:
            if self.check_tagdata(complex_mol):
                idtag = self.idtag
                system = self.system
                structure = self.structure
                positions = structure.positions
                topology = structure.topology
# begin bayly prepMD section
            #if mol.HasData(oechem.OEGetTag( self.args.PLMaskSDTag)):
            #    PLmask = json.loads(mol.GetStringData( self.args.PLMaskSDTag))
            PLmask = json.loads(complex_mol.GetStringData( self.args.PLMaskSDTag))
            # Hardwiring a strong restraint weight of 5.0 kcal/mol/ang^2 for now
            restraintWt = 5.0
            minState = plmd.RestrMin( topology, system, positions, PLmask,
                                      restraintWt, self.args.steps)
# end   bayly prepMD section
            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            complex_mol.AddData(oechem.OEGetTag('state'), output.encode(minState))
            # set up output directory if it does not already exist
            #if not os.path.exists('./output'):
            #    os.makedirs('./output')
            self.outfname = 'output/{}-minimized'.format(idtag)
            with open('output/{}-minimized.pdb'.format(idtag), 'w') as minout:
                app.PDBFile.writeFile( topology, minState.getPositions(), minout)
            self.success.emit(complex_mol)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                complex_mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(complex_mol)

class OpenMMwarmupNVTCube(OEMolComputeCube):
    title = 'Warming with strong restraints'
    description = """
    Warm up the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and warm it to the target temperature
    (default 300K) leaving all solvent free but with strong restraints
    on the protein and ligand.

    Input parameters:
    ActSiteResNumSDTag (string): the SD Tag for the whitespace-delimited protein
    active site residue numbers (integers). These residue numbers must correspond
    to residue numbers in the streamed complex.oeb.gz file.
    """
    classification = ['MDWarmup']
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    temperature = parameter.DecimalParameter(
        'temperature',
        default= 300,
        help_text="Temperature (Kelvin)")

    picoSec = parameter.DecimalParameter(
        'picoSec',
        default= 10,
        help_text="Number of picoseconds of MD")

    RestrWt = parameter.DecimalParameter(
        'RestrWt',
        default= 2,
        help_text="Restraint weight in kcal/mol per AngstromSquared")

    PLMaskSDTag = parameter.StringParameter(
        'PLMaskSDTag',
        default='OpenMM_PLmaskDict_json',
        help_text="SD Tag for the Protein-Ligand Mask for restraints")

    steps = parameter.IntegerParameter(
        'steps',
        default=1000,
        help_text="Number of MD steps")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=1000,
        help_text="Step interval for reporting data."
    )

    def check_tagdata(self, mol):
        if 'idtag' not in mol.GetData().keys():
            raise RuntimeError('Could not find idtag for molecule')
        else:
            idtag =  mol.GetData(oechem.OEGetTag('idtag'))
            if not os.path.exists('./output'):
                os.makedirs('./output')
            self.outfname = 'output/{}-simulation'.format(idtag)
            self.idtag = idtag
        if 'system' not in mol.GetData().keys():
            raise RuntimeError("Could not find system for molecule")
        else:
            sys_in = OpenMMSystemInput('sys_in')
            sys_tag = oechem.OEGetTag('system')
            system = sys_in.decode(mol.GetData(sys_tag))
            self.system = system
        if 'structure' not in mol.GetData().keys():
            raise RuntimeError('Could not find structure for molecule')
        else:
            struct_in = ParmEdStructureInput('struct_in')
            struct_tag = oechem.OEGetTag('structure')
            structure = struct_in.decode(mol.GetData(struct_tag))
            self.structure = structure

        # Check if mol has State data attached
        if 'state' in mol.GetData().keys():
            self.log.info('Found a saved State, restarting simulation')
            mol.GetData(oechem.OEGetTag('state'))
            serialized_state = mol.GetData(oechem.OEGetTag('state'))
            state = openmm.XmlSerializer.deserialize( serialized_state )
            self.state = state
            self.outfname = 'output/{}-restart'.format(self.idtag)
        else:
            self.state = None

        if not any([self.idtag, self.system, self.structure]):
            raise RuntimeError('Missing tagged generic data')
        else:
            return True

    def setReporters(self):
        from sys import stdout
        progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                            reportInterval=self.args.reporter_interval,
                                            totalSteps=self.args.steps,
                                            time=True, speed=True, progress=True,
                                            elapsedTime=True, remainingTime=True)

        state_reporter = app.StateDataReporter(self.outfname+'.log', separator="\t",
                                            reportInterval=self.args.reporter_interval,
                                            step=True,
                                            potentialEnergy=True, totalEnergy=True,
                                            volume=True, temperature=True)
        chk_reporter = app.checkpointreporter.CheckpointReporter(self.outfname+'.chk',
                                                                self.args.reporter_interval)
        import mdtraj
        traj_reporter = mdtraj.reporters.HDF5Reporter(self.outfname+'.h5', self.args.reporter_interval)
        #dcd_reporter = app.dcdreporter.DCDReporter(self.outfname+'.dcd', self.args.reporter_interval)
        self.reporters = [progress_reporter, state_reporter, traj_reporter, chk_reporter] #,dcd_reporter]
        return self.reporters

    def begin(self):
        pass

    def process(self, complex_mol, port):
        try:
            if self.check_tagdata(complex_mol):
                idtag = self.idtag
                outfname = self.outfname
                system = self.system
                structure = self.structure
                positions = structure.positions
                topology = structure.topology
            # Initialize Simulation
            integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
            simulation = app.Simulation(topology, system, integrator)
            platform = simulation.context.getPlatform().getName()
            self.log.info('Running OpenMMSimulation on Platform {}'.format(platform))

            # Check if mol has State data attached
            if self.state:
                simulation.context.setState(self.state)
            else:
                # Set initial positions and velocities then minimize
                simulation.context.setPositions(positions)
                simulation.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)
                init = simulation.context.getState(getEnergy=True)
                self.log.info('Initial energy is {}'.format(init.getPotentialEnergy()))
                self.log.info('Minimizing {} system...'.format(idtag))
                simulation.minimizeEnergy()
                st = simulation.context.getState(getPositions=True,getEnergy=True)
                self.log.info('Minimized energy is {}'.format(st.getPotentialEnergy()))
                with open('output/{}-minimized.pdb'.format(idtag), 'w') as minout:
                    app.PDBFile.writeFile(simulation.topology, st.getPositions(), minout)

            #Append Reporters to simulation
            reporters = self.setReporters()
            for rep in reporters:
                simulation.reporters.append(rep)

            self.log.info('Running {} MD steps at {}K'.format(self.args.steps, self.args.temperature))
            simulation.step(self.args.steps)
            outlog = open(outfname+'.log', 'r')
            self.log.info(outlog.read())

            # Save serialized State object
            state = simulation.context.getState(getPositions=True,
                                              getVelocities=True,
                                              getParameters=True)

            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            struct_out = ParmEdStructureOutput('struct_out')
            complex_mol.SetData(oechem.OEGetTag('system'), output.encode(system))
            complex_mol.SetData(oechem.OEGetTag('structure'), struct_out.encode(structure))
            complex_mol.AddData(oechem.OEGetTag('state'), output.encode(state))
            complex_mol.AddData(oechem.OEGetTag('log'), outlog.read())
            self.success.emit(complex_mol)
            outlog.close()

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                complex_mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(complex_mol)
