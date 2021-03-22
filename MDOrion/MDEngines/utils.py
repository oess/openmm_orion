# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
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


from orionclient.session import in_orion

from abc import ABC, abstractmethod

import parmed

import simtk

import time

import fcntl

import os

from simtk import unit

import itertools

md_keys_converter = {'OpenMM':

                         {'constraints':
                              {'None': 'None', 'Bonds2H': 'HBonds', 'Angles2H': 'HAngles', 'All-Bonds': 'AllBonds'}
                          },

                     'Gromacs':
                         {'constraints':
                              {'None': 'none', 'Bonds2H': 'h-bonds', 'Angles2H': 'h-angles', 'All-Bonds': 'all-bonds'}
                          }
                     }


class MDState(object):

    def __init__(self, parmed_structure):

        if not parmed_structure.positions:
            raise RuntimeError('Atom positions are not defined')
        else:
            # The returned object is an OpenMM Quantity with units
            self.__positions__ = parmed_structure.positions

        if parmed_structure.velocities is None:
            self.__velocities__ = None
        else:
            # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
            self.__velocities__ = parmed_structure.velocities * simtk.unit.angstrom/simtk.unit.picosecond
            # The returned object is an OpenMM Quantity with units

        if parmed_structure.box_vectors is None:
            self.__box_vectors__ = None
        else:
            self.__box_vectors__ = parmed_structure.box_vectors

    def get_positions(self):
        return self.__positions__

    def get_oe_positions(self):
        pos = self.__positions__.in_units_of(unit.angstrom) / unit.angstrom
        pos = list(itertools.chain.from_iterable(pos))
        return pos

    def get_velocities(self):
        return self.__velocities__

    def get_box_vectors(self):
        return self.__box_vectors__

    def set_positions(self, positions):
        if isinstance(positions, simtk.unit.quantity.Quantity):
            self.__positions__ = positions
        else:
            raise ValueError("It was not possible to set the positions")

        return

    def set_velocities(self, velocities):
        if isinstance(velocities, simtk.unit.quantity.Quantity):
            self.__velocities__ = velocities
        else:
            raise ValueError("It was not possible to set the velocities")

        return

    def set_box_vectors(self, box_vectors):
        if isinstance(box_vectors, simtk.unit.quantity.Quantity):
            self.__box_vectors__ = box_vectors
        else:
            raise ValueError("It was not possible to set the box vectors")

        return


class MDSimulations(ABC):
    @abstractmethod
    def __init__(self, mdstate, ff_parameters, opt):

        if not isinstance(mdstate, MDState):
            raise ValueError("{} is not a MDState Object".format(type(mdstate)))

        if not isinstance(ff_parameters, parmed.Structure):
            raise ValueError("{} is not a Parmed Structure Object".format(type(mdstate)))

        if not isinstance(opt, dict):
            raise ValueError("{} is not a dictionary".format(type(opt)))

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def update_state(self):
        pass

    @abstractmethod
    def clean_up(self):
        pass


def local_cluster(sim):

    def wrapper(*args):

        mdstate = args[0]
        ff_parameters = args[1]
        opt = args[2]

        if 'OE_VISIBLE_DEVICES' in os.environ and not in_orion():

            gpus_available_indexes = os.environ["OE_VISIBLE_DEVICES"].split(',')

            opt['Logger'].info("OE LOCAL FLOE CLUSTER OPTION IN USE")

            if 'OE_MAX' in os.environ:
                opt['OE_MAX'] = int(os.environ["OE_MAX"])
            else:
                opt['OE_MAX'] = 1

            opt['Logger'].info("OE MAX = {}".format(opt['OE_MAX']))

            while True:

                for gpu_id in gpus_available_indexes:

                    for p in range(0, opt['OE_MAX']):

                        fn = str(gpu_id) + '_' + str(p) + '.txt'

                        try:
                            with open(fn, 'a') as file:

                                fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                                # opt['Logger'].warn("LOCKED GPU ID = {} - MOL ID = {}".format(gpu_id, opt['system_id']))
                                file.write(
                                        "MD - name = {} MOL_ID = {} GPU_IDS = {} GPU_ID = {}\n".format(opt['system_title'],
                                                                                                       opt['system_id'],
                                                                                                       gpus_available_indexes,
                                                                                                       str(gpu_id)))
                                opt['gpu_id'] = str(gpu_id)

                                new_mdstate = sim(mdstate, ff_parameters, opt)

                                time.sleep(5.0)
                                # opt['Logger'].warn("UNLOCKING GPU ID = {} - MOL ID = {}".format(gpu_id, opt['system_id']))
                                fcntl.flock(file, fcntl.LOCK_UN)
                                return new_mdstate

                        except BlockingIOError:
                            time.sleep(0.1)

                        except Exception as e:  # If the simulation fails for other reasons
                            try:
                                time.sleep(5.0)
                                fcntl.flock(file, fcntl.LOCK_UN)
                            except Exception as e:
                                pass
                            raise ValueError("{} Simulation Failed".format(str(e)))
        else:
            new_mdstate = sim(*args)
            return new_mdstate

    return wrapper


@local_cluster
def md_simulation(mdstate, ff_parameters, opt):

    if opt['md_engine'] == 'OpenMM':

        from MDOrion.MDEngines.OpenMMCubes.simtools import OpenMMSimulations

        MDSim = OpenMMSimulations(mdstate, ff_parameters, opt)

        MDSim.run()

        new_mdstate = MDSim.update_state()

        MDSim.clean_up()

        return new_mdstate

    if opt['md_engine'] == 'Gromacs':

        from MDOrion.MDEngines.Gromacs.simtools import GromacsSimulations

        MDSim = GromacsSimulations(mdstate, ff_parameters, opt)

        MDSim.run()

        new_mdstate = MDSim.update_state()

        MDSim.clean_up()

        return new_mdstate

    else:
        raise ValueError("The selected MD engine is not currently supported: {}".format(opt['md_engine']))
