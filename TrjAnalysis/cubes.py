from floe.api import (parameter,
                      ParallelMixin)

from cuberecord import (OERecordComputeCube,
                        OEField)
from cuberecord.constants import DEFAULT_MOL_NAME
from datarecord import (Types,
                        Meta,
                        ColumnMeta)

import traceback
from OpenMMCubes.utils import ParmedData
from oeommtools.utils import split
from TrjAnalysis.utils import BoundingBox
import sstmap as sm

from openeye import oechem
import mdtraj as md
import numpy as np
from openeye.oegrid import OEScalarGrid, OEWriteGrid

from tempfile import TemporaryDirectory
import os


class SSTMapGistCube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "Water Thermodynamics by Using Gist in SSTMap"
    description = """
    SSTMap performs Water Thermodynamics analysis.
    SSTMaps supports hydration site analysis (HSA) 
    and Grid Inhomogeneous Solvation Theory (GIST). 
    
    SSTMap has been developed at Kurtzman Lab Lehman College 
    For more details, please visit 
    sstmap.org @ https://github.com/KurtzmanLab/SSTMap
    """
    classification = ["SSTMap Analysis"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    center_grid_ligand = parameter.BooleanParameter(
        'center_grid_ligand',
        default=True,
        help_text='Center the grid around the ligand instead of the protein'

    )

    trj_fn = parameter.StringParameter(
        'trj_fn',
        default='trj.dcd',
        help_text='Trajectory file name'
    )

    grid_resolution = parameter.DecimalParameter(
        'grid_resolution',
        default=0.5,
        help_text='Grid resolution in A'
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            field_system = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                   meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))

            if not record.has_value(field_system):
                self.log.warn("Missing molecule '{}' field".format(field_system.get_name()))
                self.failure.emit(record)
                return

            system = record.get_value(field_system)

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = split(system, ligand_res_name='LIG')

            self.log.info("\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))

            if opt['center_grid_ligand'] and ligand.NumAtoms():
                bb = BoundingBox(ligand, save_pdb=False, scale_factor=1.4)
            elif protein.NumAtoms():
                bb = BoundingBox(protein, save_pdb=False, scale_factor=1.4)
            else:
                oechem.OEThrow.Error("Protein and Ligand molecules have not been detected")

            field_parmed = OEField("Parmed", ParmedData)

            if not record.has_value(field_parmed):
                self.log.warn("Missing molecule '{}' field".format(field_parmed.get_name()))
                self.failure.emit(record)

            parmed_structure = record.get_value(field_parmed)

            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Temporary Output Directory {}".format(output_directory))

                system_fn = os.path.join(output_directory, "system.gro")

                system_top_fn = system_fn.split(".")[0]+'.top'

                parmed_structure.save(system_fn, overwrite=True)
                parmed_structure.save(system_top_fn, overwrite=True)

                # Grid Center
                center = (bb[0]+bb[1])*0.5

                trj = md.load(opt['trj_fn'], top=system_fn)

                # Number of voxels in each dimension
                dims = (np.ceil(bb[1]) - np.floor(bb[0]))/opt['grid_resolution'] + 2
                nx = int(dims[0])
                ny = int(dims[1])
                nz = int(dims[2])

                # Define the Gist Grid
                gist = sm.GridWaterAnalysis(
                            system_fn,
                            opt['trj_fn'],
                            start_frame=0, num_frames=trj.n_frames,
                            grid_center=center,
                            grid_dimensions=[nx,
                                             ny,
                                             nz],
                            grid_resolution=[opt['grid_resolution'],
                                             opt['grid_resolution'],
                                             opt['grid_resolution']],
                            prefix="system", supporting_file=system_top_fn)

                # Info Gist
                gist.print_system_summary()

                # Run Gist
                gist.calculate_grid_quantities()

                gist.print_calcs_summary()
                # gist.write_data()

                # Write out .dx files
                # gist.generate_dx_files(prefix=output_directory)

                # Extract info from gist
                quantities = ['index', 'x', 'y', 'z',
                              'N_wat', 'g_O', 'g_H',
                              'TS_tr_dens', 'TS_tr_norm',
                              'TS_or_dens', 'TS_or_norm',
                              'dTSsix-dens', 'dTSsix-norm',
                              'E_sw_dens', 'E_sw_norm', 'E_ww_dens', 'Eww_norm',
                              'E_ww_nbr_dens', 'E_ww_nbr_norm',
                              'N_nbr_dens', 'N_nbr_norm',
                              'f_hb_dens', 'f_hb_norm',
                              'N_hb_sw_dens', 'N_hb_sw_norm', 'N_hb_ww_dens', 'N_hb_ww_norm',
                              'N_don_sw_dens', 'N_don_sw_norm', 'N_acc_sw_dens', 'N_acc_sw_norm',
                              'N_don_ww_dens', 'N_don_ww_norm', 'N_acc_ww_dens', 'N_acc_ww_norm']

                indices = range(0, 35)
                index_quantity_map = dict(zip(quantities, indices))

                # Extract voxel coordinates
                x_gist_voxels = gist.voxeldata[:, index_quantity_map["x"]]
                y_gist_voxels = gist.voxeldata[:, index_quantity_map["y"]]
                z_gist_voxels = gist.voxeldata[:, index_quantity_map["z"]]

                # Extract Structural and thermodynamics water info
                grid_gO = gist.voxeldata[:, index_quantity_map["g_O"]]
                grid_E_sw_dens = gist.voxeldata[:, index_quantity_map["E_sw_dens"]]
                grid_E_ww_dens = gist.voxeldata[:, index_quantity_map["E_ww_dens"]]

                # Generate the OEGrids
                oegrid_gO = OEScalarGrid(nx, ny, nz,
                                         center[0],
                                         center[1],
                                         center[2],
                                         opt['grid_resolution'])

                oegrid_Esw_dens = OEScalarGrid(nx, ny, nz,
                                               center[0],
                                               center[1],
                                               center[2],
                                               opt['grid_resolution'])

                oegrid_Eww_dens = OEScalarGrid(nx, ny, nz,
                                               center[0],
                                               center[1],
                                               center[2],
                                               opt['grid_resolution'])

                # Set OEGrid Values from Gist voxels
                for i in range(oegrid_gO.GetSize()):
                    x = x_gist_voxels[i]
                    y = y_gist_voxels[i]
                    z = z_gist_voxels[i]
                    oegrid_gO.SetValue(x, y, z, grid_gO[i])
                    oegrid_Esw_dens.SetValue(x, y, z, grid_E_sw_dens[i])
                    oegrid_Eww_dens.SetValue(x, y, z, grid_E_ww_dens[i])

                # Write out OEGrids
                # OEWriteGrid("g_O.grd", oegrid_gO)
                # OEWriteGrid("E_sw_dens.grd", oegrid_Esw_dens)
                # OEWriteGrid("E_ww_dens.grd", oegrid_Eww_dens)

                field_gO = OEField("gO", Types.Chem.Grid)
                field_esw = OEField("Esw", Types.Chem.Grid)
                field_eww = OEField("Eww", Types.Chem.Grid)

                record.set_value(field_gO, oegrid_gO)
                record.set_value(field_esw, oegrid_Esw_dens)
                record.set_value(field_eww, oegrid_Eww_dens)

            # Emit the ligand
            self.success.emit(record)

        except Exception as e:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class SSTMapHSACube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "Water Thermodynamics by Using HSA in SSTMap"
    description = """
    SSTMap performs Water Thermodynamics analysis.
    SSTMaps supports hydration site analysis (HSA) 
    and Grid Inhomogeneous Solvation Theory (GIST). 

    SSTMap has been developed at Kurtzman Lab Lehman College 
    For more details, please visit 
    sstmap.org @ https://github.com/KurtzmanLab/SSTMap
    """
    classification = ["SSTMap Analysis"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    trj_fn = parameter.StringParameter(
        'trj_fn',
        default='trj.dcd',
        help_text='Trajectory file name'
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            field_system = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                   meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))

            if not record.has_value(field_system):
                self.log.warn("Missing molecule '{}' field".format(field_system.get_name()))
                self.failure.emit(record)
                return

            system = record.get_value(field_system)

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = split(system, ligand_res_name='LIG')

            self.log.info("\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))
            field_parmed = OEField("Parmed", ParmedData)

            if not record.has_value(field_parmed):
                self.log.warn("Missing molecule '{}' field".format(field_parmed.get_name()))
                self.failure.emit(record)

            parmed_structure = record.get_value(field_parmed)

            # with TemporaryDirectory() as output_directory:

                #opt['Logger'].info("Temporary Output Directory {}".format(output_directory))

                # system_fn = os.path.join(output_directory, "system.gro")

            system_fn = "system.gro"

            system_top_fn = system_fn.split(".")[0] + '.top'

            ligand_fn = "ligand.pdb"

            with oechem.oemolostream(ligand_fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, ligand)

            parmed_structure.save(system_fn, overwrite=True)
            parmed_structure.save(system_top_fn, overwrite=True)

            trj = md.load(opt['trj_fn'], top=system_fn)

            # Set SSTMap in HSA mode
            hsa = sm.SiteWaterAnalysis(system_fn,
                                       opt["trj_fn"],
                                       start_frame=0, num_frames=trj.n_frames,
                                       supporting_file=system_top_fn,
                                       hsa_region_radius=5.0,
                                       ligand_file=ligand_fn,
                                       prefix="system")

            # Info HSA
            hsa.print_system_summary()

            # Run HSA

            # Writing starting cluster file
            hsa.initialize_hydration_sites()

            # Write hydration sites and probable water configuration
            hsa.calculate_site_quantities()

            # hsa.write_calculation_summary()
            # hsa.write_data()

        except Exception as e:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return