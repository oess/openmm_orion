import io, os, time, traceback
import numpy as np
import json
from openeye import oechem

from simtk import unit, openmm
from simtk.openmm import app

import OpenMMCubes.simtools as simtools
import OpenMMCubes.utils as utils
from OpenMMCubes import plmd

from floe.api import OEMolComputeCube, parameter, ParallelOEMolComputeCube

def ExtractOpenMMData( mol):
    if 'idtag' not in mol.GetData().keys():
        raise RuntimeError('Could not find idtag for molecule')
    else:
        idtag =  mol.GetData(oechem.OEGetTag('idtag'))
    if 'system' not in mol.GetData().keys():
        raise RuntimeError("Could not find system for molecule")
    else:
        sys_in = OpenMMSystemInput('sys_in')
        sys_tag = oechem.OEGetTag('system')
        system = sys_in.decode(mol.GetData(sys_tag))
    if 'structure' not in mol.GetData().keys():
        raise RuntimeError('Could not find structure for molecule')
    else:
        struct_in = ParmEdStructureInput('struct_in')
        struct_tag = oechem.OEGetTag('structure')
        structure = struct_in.decode(mol.GetData(struct_tag))
        positions = structure.positions
        topology = structure.topology
    # Check if mol has State data attached
    if 'state' in mol.GetData().keys():
        mol.GetData(oechem.OEGetTag('state'))
        serialized_state = mol.GetData(oechem.OEGetTag('state'))
        state = openmm.XmlSerializer.deserialize( serialized_state )
    else:
        state = None
    # Check if mol has protein-ligand mask data attached
    if 'OpenMM_PLmaskDict_json' in mol.GetData().keys():
        PLmask = json.loads( mol.GetStringData( 'OpenMM_PLmaskDict_json'))
    else:
        PLmask = None

    if not any([idtag, system, structure]):
        raise RuntimeError('Missing tagged generic data')
    else:
        OpenMMstuff = { 'idtag': idtag,
                        'system': system,
                        'topology': topology,
                        'positions': positions,
                        'state': state,
                        'PLmask': PLmask }
        return OpenMMstuff

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
    #intake = CustomMoleculeInputPort('intake')
    #success = CustomMoleculeOutputPort('success')

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
    #intake = CustomMoleculeInputPort('intake')
    #success = CustomMoleculeOutputPort('success')

    steps = parameter.IntegerParameter(
        'steps',
        default=100,
        help_text="Number of minimization steps")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=5.0,
        help_text="Restraint weight for xyz atom restraints")

    def begin(self):
        if not os.path.exists('./output'):
            os.makedirs('./output')
        return

    def process(self, complex_mol, port):
        try:
            openmmStuff = ExtractOpenMMData( complex_mol)
            argsDict = vars( self.args)
            minState = plmd.RestrMin( openmmStuff, argsDict)
            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            complex_mol.SetData(oechem.OEGetTag('state'), output.encode(minState))
            outfname = 'output/{}-minimized'.format(openmmStuff['idtag'])
            with open(outfname+'.pdb', 'w') as minout:
                app.PDBFile.writeFile( openmmStuff['topology'], minState.getPositions(), minout)
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
      picosec (int): Number of picoseconds to warm up the complex.
      temperature (decimal): target final temperature after warming.
      restraintWt (decimal): strength in kcal/mol/ang^2 for xyz atom restraints.
    """
    classification = ['MDWarmup']
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    #intake = CustomMoleculeInputPort('intake')
    #success = CustomMoleculeOutputPort('success')

    temperature = parameter.DecimalParameter(
        'temperature',
        default= 300,
        help_text="Temperature (Kelvin)")

    picosec = parameter.DecimalParameter(
        'picosec',
        default= 10,
        help_text="Number of picoseconds of MD")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight in kcal/mol/ang^2 for xyz atom restraints")

    def begin(self):
        if not os.path.exists('./output'):
            os.makedirs('./output')
        return

    def process(self, complex_mol, port):
        try:
            openmmStuff = ExtractOpenMMData( complex_mol)
            argsDict = vars( self.args)
            warmState = plmd.RestrWarmupNVT( openmmStuff, argsDict)

            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            complex_mol.SetData(oechem.OEGetTag('state'), output.encode(warmState))
            outfname = 'output/{}-warmup'.format(openmmStuff['idtag'])
            with open(outfname+'.pdb', 'w') as out:
                app.PDBFile.writeFile( openmmStuff['topology'], warmState.getPositions(), out)
            self.success.emit(complex_mol)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                complex_mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(complex_mol)

class OpenMMequilCube(OEMolComputeCube):
    title = 'Equilibratiing with restraints'
    description = """
    Equilibrate the warmed up solvated protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the warmed up solvated protein:ligand complex and equilibrate it at
    the target temperature leaving all solvent free but with restraints
    on the protein and ligand.

    Input parameters:
      picosec (int): Number of picoseconds to warm up the complex.
      temperature (decimal): target temperature for equilibration.
      restraintWt (decimal): strength in kcal/mol/ang^2 for xyz atom restraints.
      snapFreq (int): frequency (in picoseconds) for taking snapshots.
    """
    classification = ['MDWarmup']
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    #intake = CustomMoleculeInputPort('intake')
    #success = CustomMoleculeOutputPort('success')

    picosec = parameter.DecimalParameter(
        'picosec',
        default= 10,
        help_text="Number of picoseconds of MD")

    snapFreq = parameter.DecimalParameter(
        'snapFreq',
        default= 0,
        help_text="frequency (in picoseconds) for taking snapshots")

    restraintType = parameter.StringParameter(
        'restraintType',
        default='NonSolventNonH',
        choices=[ 'NonSolventNonH', 'ProteinNonH', 'ProteinCAlpha', 'LigandNonH',
                  'CAlphaLigandNonH', 'None'],
        help_text= 'Which kind of atoms get xyz restraints')

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight in kcal/mol/ang^2 for xyz atom restraints")

    temperature = parameter.DecimalParameter(
        'temperature',
        default= 300,
        help_text="Temperature (Kelvin)")

    label = parameter.StringParameter(
        'label',
        default='_equil',
        help_text= 'label to add to filenames from this cube')

    def begin(self):
        self.argsDict = vars( self.args)
        return

    def process(self, complex_mol, port):
        try:
            #openmmStuff = ExtractOpenMMData( complex_mol)
            openmmStuff = utils.PackageOEMol.unpack( complex_mol)
            self.argsDict['outfname'] = '{}'.format(openmmStuff['IDTag'])+self.args.label
            equilSimuln = plmd.RestrEquil( openmmStuff, **self.argsDict)

            # Attach openmm objects to mol, emit to output
            #output = OpenMMSystemOutput('output')
            #complex_mol.SetData(oechem.OEGetTag('state'), output.encode(equilState))
            packedmol = utils.PackageOEMol.pack(complex_mol, simulation)
            packedmol.SetData(oechem.OEGetTag( 'outfname'), self.argsDict['outfname'])
            utils.PackageOEMol.dump(
                    packedmol, outfname=self.argsDict['outfname'], tarxz=True )
            self.success.emit(complex_mol)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                complex_mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(complex_mol)
