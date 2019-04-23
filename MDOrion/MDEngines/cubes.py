# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.


import traceback

from floe.api import (ParallelMixin,
                      parameter)

from cuberecord import OERecordComputeCube

from MDOrion.Standards import (MDStageTypes,
                               MDEngines)

from MDOrion.Standards.mdrecord import MDDataRecord

from MDOrion.MDEngines.utils import md_simulation

import copy

import textwrap

import os


class MDMinimizeCube(ParallelMixin, OERecordComputeCube):
    title = 'Minimization Cube'

    version = "0.1.0"
    classification = [["MD Simulations"]]
    tags = ['OpenMM', 'Gromacs', 'Minimization']

    description = """
    This cube performs energy minimization on the provided system. The system 
    must have been parametrized by the Force Field cube and the system Parmed 
    object must be present on the input record. In addition, a system identification 
    number must be present on the input record as well. This can be accomplished 
    by using the “ID Setting Cube”. The system minimization is performed by 
    the selected MD engine, currently OpenMM and Gromacs only. Restraints 
    and constraints can be used as well. Currently implicit solvent models can 
    be used in OpenMM only. The cube requires a record as input and produces 
    a new record with the minimized system.
    
    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the systems to minimize

    Output:
    -------
    oechem.OEDataRecord - Streamed-out of records with the systems energy 
    minimized.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda9"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    steps = parameter.IntegerParameter(
        'steps',
        default=0,
        help_text="""Number of minimization steps.
                  If 0 the minimization will continue
                  until convergence""")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text=""""Mask selection to apply harmonic restraints. 
        Possible keywords are: ligand, protein, water, ions, 
        ca_protein, cofactors. The selection can be refined 
        by using logical tokens: not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=5.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD
        simulation. Possible keywords are: ligand, protein, water,
        ions, ca_protein, cofactors. The selection can be refined by
        using logical tokens: not, noh, and, or, diff, around. NOTE:
        Not currently implemented in Gromacs""")

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if the non-bonded method is NoCutoff""")

    constraints = parameter.StringParameter(
        'constraints',
        default='H-Bonds',
        choices=['None', 'H-Bonds', 'H-Angles', 'All-Bonds'],
        help_text="""None, H-Bonds, H-Angles, or All-Bonds
        Which type of constraints to add to the system.
        None means no bonds are constrained.
        H-Bonds means bonds with hydrogen are constrained, etc.""")

    implicit_solvent = parameter.StringParameter(
        'implicit_solvent',
        default='None',
        choices=['None', 'HCT', 'OBC1', 'OBC2', 'GBn', 'GBn2'],
        help_text="Implicit Solvent Model. NOTE:"
                  "Not currently implemented in Gromacs")

    center = parameter.BooleanParameter(
        'center',
        default=True,
        description='Center the system to the OpenMM and Gromacs unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity. NOTE:'
                    'Not currently implemented in Gromacs')

    suffix = parameter.StringParameter(
        'suffix',
        default='min',
        help_text='Filename suffix for output simulation files')

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        description='On enables Hydrogen Mass Repartitioning. NOTE:'
                    'Not currently implemented in Gromacs')

    save_md_stage = parameter.BooleanParameter(
        'save_md_stage',
        default=True,
        help_text="""Save the MD simulation stage. If True the MD,
           simulation data will be appended to the md simulation stages 
           otherwise the last MD stage will be overwritten""")

    md_engine = parameter.StringParameter(
        'md_engine',
        default='OpenMM',
        choices=['OpenMM', 'Gromacs'],
        help_text='Select the MD available engine')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'min'
        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt['CubeTitle'] = self.title

            # Logger string
            str_logger = '-'*32 + ' MIN CUBE PARAMETERS ' + '-'*32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt['CubeTitle'])

            for k, v in sorted(self.parameters().items()):
                tmp_default = copy.deepcopy(v)

                if v.default is None:
                    tmp_default.default = 'None'
                elif isinstance(v, parameter.BooleanParameter):
                    if v.default:
                        tmp_default.default = 'True'
                    else:
                        tmp_default.default = 'False'
                else:
                    tmp_description = textwrap.fill(" ".join(v.description.split()),
                                                    subsequent_indent=' ' * 39, width=80)
                    str_logger += "\n{:<25} = {:<10} {}".format(k,
                                                                getattr(self.args, tmp_default.name),
                                                                tmp_description)

            str_logger += "\n{:<25} = {:<10}".format("Simulation Type", opt['SimType'])

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            system = mdrecord.get_primary

            if not mdrecord.has_title:
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = mdrecord.get_title

            opt['system_title'] = system_title
            opt['system_id'] = mdrecord.get_id

            system = mdrecord.get_stage_topology()
            mdstate = mdrecord.get_stage_state()

            opt['out_directory'] = mdrecord.cwd
            opt['molecule'] = system
            opt['str_logger'] = str_logger
            opt['Logger'].info('[{}] MINIMIZING System: {}'.format(opt['CubeTitle'], system_title))

            # Extract the Parmed structure and synchronize it with the last MD stage state
            parmed_structure = mdrecord.get_parmed(sync_stage_name='last')

            # Run the MD simulation
            new_mdstate = md_simulation(mdstate, parmed_structure, opt)

            # Update the system coordinates
            system.SetCoords(new_mdstate.get_oe_positions())
            mdrecord.set_primary(system)

            data_fn = os.path.basename(mdrecord.cwd) + '_' + opt['system_title'] + '_' + str(opt['system_id']) + '-' + opt['suffix'] + '.tar.gz'

            if not mdrecord.add_new_stage(self.title,
                                          MDStageTypes.MINIMIZATION,
                                          system,
                                          new_mdstate,
                                          data_fn,
                                          append=opt['save_md_stage'],
                                          log=opt['str_logger']):

                raise ValueError("Problems adding the new Minimization Stage")

            # Synchronize the added Parmed structure with the last MD stage state
            # mdrecord.set_parmed(parmed_structure, sync_stage_name='last')

            self.success.emit(mdrecord.get_record)

            del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class MDNvtCube(ParallelMixin, OERecordComputeCube):
    title = 'NVT Cube'
    version = "0.1.0"
    classification = [["MD Simulations"]]
    tags = ['Gromacs', 'OpenMM', 'NVT']

    description = """
    This cube performs MD simulation in the NVT ensemble on the provided system. 
    The system must have been parametrized by the Force Field cube and the system Parmed 
    object must be present on the input record. In addition, a system identification 
    number must be present on the input record as well. This can be accomplished 
    by using the “ID Setting Cube”. The NVT MD simulation is performed by the selected 
    MD engine, currently OpenMM and Gromacs only. Restraints and constraints can be 
    used as well. Currently implicit solvent models can be used in OpenMM only. 
    The cube requires a record as input and produces a new record with the time evolved 
    system. The total sampling time can be set by using the “time” cube parameter 
    while the trajectory snapshots can be set by using the “trajectory_interval” cube 
    parameter.


    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the systems to NVT sample

    Output:
    -------
    oechem.OEDataRecord - Streamed-out of records with the systems NVT time evolved
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda9"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    time = parameter.DecimalParameter(
        'time',
        default=0.01,
        help_text="NVT simulation time in nanoseconds")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text=""""Mask selection to apply harmonic restraints. 
        Possible keywords are: ligand, protein, water, ions, 
        ca_protein, cofactors. The selection can be refined 
        by using logical tokens: not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD
        simulation. Possible keywords are: ligand, protein, water,
        ions, ca_protein, cofactors. The selection can be refined by
        using logical tokens: not, noh, and, or, diff, around. NOTE:
        Not currently implemented in Gromacs""")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff.""")

    constraints = parameter.StringParameter(
        'constraints',
        default='H-Bonds',
        choices=['None', 'H-Bonds', 'H-Angles', 'All-Bonds'],
        help_text="""None, H-Bonds, H-Angles, or All-Bonds
        Which type of constraints to add to the system.
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained, etc.""")

    implicit_solvent = parameter.StringParameter(
        'implicit_solvent',
        default='None',
        choices=['None', 'HCT', 'OBC1', 'OBC2', 'GBn', 'GBn2'],
        help_text="Implicit Solvent Model. NOTE:"
                  "Not currently implemented in Gromacs")

    trajectory_interval = parameter.DecimalParameter(
        'trajectory_interval',
        default=0.0,
        help_text="""Time interval for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""")

    reporter_interval = parameter.DecimalParameter(
        'reporter_interval',
        default=0.0,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""")

    suffix = parameter.StringParameter(
        'suffix',
        default='nvt',
        help_text='Filename suffix for output simulation files')

    center = parameter.BooleanParameter(
        'center',
        default=False,
        help_text='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        help_text='Increase log file verbosity. NOTE: '
                  'Not currently implemented in Gromacs')

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        help_text='On enables Hydrogen Mass Repartitioning. NOTE: '
                  'Not currently implemented in Gromacs')

    save_md_stage = parameter.BooleanParameter(
        'save_md_stage',
        default=False,
        help_text="""Save the MD simulation stage. If True the MD,
           simulation data will be appended to the md simulation stages 
           otherwise the last MD stage will be overwritten""")

    md_engine = parameter.StringParameter(
        'md_engine',
        default='OpenMM',
        choices=['OpenMM', 'Gromacs'],
        help_text='Select the MD available engine')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'nvt'

        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt['CubeTitle'] = self.title

            # Logger string
            str_logger = '-'*32 + ' NVT CUBE PARAMETERS ' + '-'*32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt['CubeTitle'])

            for k, v in sorted(self.parameters().items()):
                tmp_default = copy.deepcopy(v)

                if v.default is None:
                    tmp_default.default = 'None'
                elif isinstance(v, parameter.BooleanParameter):
                    if v.default:
                        tmp_default.default = 'True'
                    else:
                        tmp_default.default = 'False'
                else:
                    tmp_description = textwrap.fill(" ".join(v.description.split()),
                                                    subsequent_indent=' ' * 39, width=80)
                    str_logger += "\n{:<25} = {:<10} {}".format(k,
                                                                getattr(self.args, tmp_default.name),
                                                                tmp_description)

            str_logger += "\n{:<25} = {:<10}".format("Simulation Type", opt['SimType'])

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            system = mdrecord.get_primary

            if not mdrecord.has_title:
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = mdrecord.get_title

            opt['system_title'] = system_title
            opt['system_id'] = mdrecord.get_id

            system = mdrecord.get_stage_topology()
            mdstate = mdrecord.get_stage_state()

            opt['out_directory'] = mdrecord.cwd
            opt['molecule'] = system
            opt['str_logger'] = str_logger
            opt['Logger'].info('[{}] START NVT SIMULATION: {}'.format(opt['CubeTitle'], system_title))

            opt['out_fn'] = os.path.basename(opt['out_directory']) + '_' + \
                            opt['system_title'] + '_' + \
                            str(opt['system_id']) + '-' + \
                            opt['suffix']

            # Trajectory file name if any generated
            opt['trj_fn'] = opt['out_fn'] + '_' + 'traj.tar.gz'

            # Extract the Parmed structure and synchronize it with the last MD stage state
            parmed_structure = mdrecord.get_parmed(sync_stage_name='last')

            # Run the MD simulation
            new_mdstate = md_simulation(mdstate, parmed_structure, opt)

            # Update the system coordinates
            system.SetCoords(new_mdstate.get_oe_positions())
            mdrecord.set_primary(system)

            # Trajectory
            if opt['trajectory_interval']:
                trajectory_fn = opt['trj_fn']
                if opt['md_engine'] == MDEngines.OpenMM:
                    trajectory_engine = MDEngines.OpenMM
                else:
                    trajectory_engine = MDEngines.Gromacs
            else:  # Empty Trajectory
                trajectory_fn = None
                trajectory_engine = None

            data_fn = opt['out_fn']+'.tar.gz'

            if not mdrecord.add_new_stage(self.title,
                                          MDStageTypes.NVT,
                                          system,
                                          new_mdstate,
                                          data_fn,
                                          append=opt['save_md_stage'],
                                          log=opt['str_logger'],
                                          trajectory_fn=trajectory_fn,
                                          trajectory_engine=trajectory_engine,
                                          trajectory_orion_ui=opt['system_title'] + '_' + str(opt['system_id']) + '-' + opt['suffix']+'.tar.gz'
                                          ):

                raise ValueError("Problems adding in the new NVT Stage")

            # Synchronize the added Parmed structure with the last MD stage state
            # mdrecord.set_parmed(parmed_structure, sync_stage_name='last')

            self.success.emit(mdrecord.get_record)

            del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class MDNptCube(ParallelMixin, OERecordComputeCube):
    title = 'NPT Cube'
    version = "0.1.0"
    classification = [['MD Simulations']]
    tags = ['Gromacs', 'OpenMM', 'NPT']

    description = """
    This cube performs MD simulation in the NPT ensemble on the provided system. 
    The system must have been parametrized by the Force Field cube and the system Parmed 
    object must be present on the input record. In addition, a system identification 
    number must be present on the input record as well. This can be accomplished 
    by using the “ID Setting Cube”. The NPT MD simulation is performed by the selected 
    MD engine, currently OpenMM and Gromacs only. Restraints and constraints can be 
    used as well. Currently implicit solvent models can be used in OpenMM only. 
    The cube requires a record as input and produces a new record with the time evolved 
    system. The total sampling time can be set by using the “time” cube parameter 
    while the trajectory snapshots can be set by using the “trajectory_interval” cube 
    parameter.


    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the systems to NVT sample

    Output:
    -------
    oechem.OEDataRecord - Streamed-out of records with the systems NPT time evolved
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda9"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    time = parameter.DecimalParameter(
        'time',
        default=0.01,
        help_text="NPT simulation time in nanoseconds")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text=""""Mask selection to apply harmonic restraints. 
        Possible keywords are: ligand, protein, water, ions, 
        ca_protein, cofactors. The selection can be refined 
        by using logical tokens: not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD simulation.
        Possible keywords are: ligand, protein, water, ions, ca_protein,
        cofactors. The selection can be refined by using logical tokens:
        not, noh, and, or, diff, around. Not currently implemented in Gromacs""")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff""")

    constraints = parameter.StringParameter(
        'constraints',
        default='H-Bonds',
        choices=['None', 'H-Bonds', 'H-Angles', 'All-Bonds'],
        help_text="""None, H-Bonds, H-Angles, or All-Bonds
        Which type of constraints to add to the system.
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained, etc.""")

    implicit_solvent = parameter.StringParameter(
        'implicit_solvent',
        default='None',
        choices=['None', 'HCT', 'OBC1', 'OBC2', 'GBn', 'GBn2'],
        help_text="Implicit Solvent Model. Not Currently implemented in Gromacs")

    trajectory_interval = parameter.DecimalParameter(
        'trajectory_interval',
        default=0.0,
        help_text="""Time interval for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""")

    reporter_interval = parameter.DecimalParameter(
        'reporter_interval',
        default=0.0,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""")

    suffix = parameter.StringParameter(
        'suffix',
        default='npt',
        help_text='Filename suffix for output simulation files')

    center = parameter.BooleanParameter(
        'center',
        default=False,
        help_text='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        help_text='Increase log file verbosity')

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        help_text='On enables Hydrogen Mass Repartitioning. Not currently implemented in Gromacs')

    save_md_stage = parameter.BooleanParameter(
        'save_md_stage',
        default=False,
        help_text="""Save the MD simulation stage. If True the MD,
           simulation data will be appended to the md simulation stages 
           otherwise the last MD stage will be overwritten""")

    md_engine = parameter.StringParameter(
        'md_engine',
        default='OpenMM',
        choices=['OpenMM', 'Gromacs'],
        help_text='Select the MD available engine')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'npt'

        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt['CubeTitle'] = self.title
            # Logger string
            str_logger = '-'*32 + ' NPT CUBE PARAMETERS ' + '-'*32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt['CubeTitle'])

            for k, v in sorted(self.parameters().items()):
                tmp_default = copy.deepcopy(v)

                if v.default is None:
                    tmp_default.default = 'None'
                elif isinstance(v, parameter.BooleanParameter):
                    if v.default:
                        tmp_default.default = 'True'
                    else:
                        tmp_default.default = 'False'
                else:
                    tmp_description = textwrap.fill(" ".join(v.description.split()),
                                                    subsequent_indent=' ' * 39, width=80)
                    str_logger += "\n{:<25} = {:<10} {}".format(k,
                                                                getattr(self.args, tmp_default.name),
                                                                tmp_description)

            str_logger += "\n{:<25} = {:<10}".format("Simulation Type", opt['SimType'])

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            system = mdrecord.get_primary

            if not mdrecord.has_title:
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = mdrecord.get_title

            opt['system_title'] = system_title
            opt['system_id'] = mdrecord.get_id

            system = mdrecord.get_stage_topology()
            mdstate = mdrecord.get_stage_state()

            opt['out_directory'] = mdrecord.cwd
            opt['molecule'] = system
            opt['str_logger'] = str_logger
            opt['Logger'].info('[{}] START NPT SIMULATION: {}'.format(opt['CubeTitle'], system_title))

            opt['out_fn'] = os.path.basename(opt['out_directory']) + '_' + \
                            opt['system_title'] + '_' + \
                            str(opt['system_id']) + '-' + \
                            opt['suffix']

            # Trajectory file name if any generated
            opt['trj_fn'] = opt['out_fn'] + '_' + 'traj.tar.gz'

            # Extract the Parmed structure and synchronize it with the last MD stage state
            parmed_structure = mdrecord.get_parmed(sync_stage_name='last')

            # Run the MD simulation
            new_mdstate = md_simulation(mdstate, parmed_structure, opt)

            # Update the system coordinates
            system.SetCoords(new_mdstate.get_oe_positions())
            mdrecord.set_primary(system)

            # Trajectory
            if opt['trajectory_interval']:
                trajectory_fn = opt['trj_fn']
                if opt['md_engine'] == MDEngines.OpenMM:
                    trajectory_engine = MDEngines.OpenMM
                else:
                    trajectory_engine = MDEngines.Gromacs

            else:  # Empty Trajectory
                trajectory_fn = None
                trajectory_engine = None

            data_fn = opt['out_fn'] + '.tar.gz'

            if not mdrecord.add_new_stage(self.title,
                                          MDStageTypes.NPT,
                                          system,
                                          new_mdstate,
                                          data_fn,
                                          append=opt['save_md_stage'],
                                          log=opt['str_logger'],
                                          trajectory_fn=trajectory_fn,
                                          trajectory_engine=trajectory_engine,
                                          trajectory_orion_ui=opt['system_title'] + '_' + str(opt['system_id']) + '-' + opt['suffix']+'.tar.gz'
                                          ):

                raise ValueError("Problems adding in the new NPT Stage")

            # Synchronize the added Parmed structure with the last MD stage state
            # mdrecord.set_parmed(parmed_structure, sync_stage_name='last')

            self.success.emit(mdrecord.get_record)

            del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return
