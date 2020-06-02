from orionplatform.mixins import RecordPortsMixin

from floe.api import (ParallelMixin,
                      parameters,
                      ComputeCube)

from MDOrion.Standards import Fields, MDStageNames

from oeommtools import utils as oeutils

from floereport import FloeReport, LocalFloeReport

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

from os import environ

import MDOrion.TrjAnalysis.utils as utl

import MDOrion.TrjAnalysis.TrajMMPBSA_utils as mmpbsa

from MDOrion.TrjAnalysis.water_utils import nmax_waters

import oetrajanalysis.OETrajBasicAnalysis_utils as oetrjutl

import ensemble2img

from tempfile import TemporaryDirectory

from openeye import oechem

import oetrajanalysis.Clustering_utils as clusutl

from openeye import oedepict

import os

import traceback

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from MDOrion.Standards.mdrecord import MDDataRecord

from MDOrion.TrjAnalysis.TrajAnFloeReport_utils import (_clus_floe_report_header,
                                                        _clus_floe_report_header2,
                                                        _clus_floe_report_midHtml0,
                                                        _clus_floe_report_midHtml1,
                                                        _clus_floe_report_midHtml2,
                                                        _clus_floe_report_stripPlots,
                                                        _clus_floe_report_Trailer,
                                                        trim_svg,
                                                        MakeClusterInfoText)

# import traceback
#
# from floe.api import ParallelMixin, parameter
#
# from cuberecord import OERecordComputeCube
#
# from cuberecord.ports import RecordInputPort
#
# from datarecord import (OEField,
#                         Types,
#                         OERecord)
#
#
# from tempfile import TemporaryDirectory
#
# import os
#
# from oeommtools.utils import split
#
# from Standards import Fields, MDStageTypes
#
# import mdtraj as md
#
# import sstmap as sm
#
# from MDEngines import utils as omm_utils
#
# from floe.constants import *
#
# from openeye import oechem, oegrid
#
# from TrjAnalysis import sstmap_utils
#
# from TrjAnalysis.sstmap_utils import GISTFields
#
# import copy as cp
#
# import shutil
#
# from orionclient.session import in_orion
#
# from shutil import copyfile
#
#
# class SSTMapHsa(ParallelMixin, OERecordComputeCube):
#
#     # version = "0.1.4"
#
#     title = "SSTMAP HSA Analysis"
#
#     description = """
#         SSTMap performs Water Thermodynamics analysis.
#         SSTMaps supports hydration site analysis (HSA)
#         and Grid Inhomogeneous Solvation Theory (GIST).
#
#         SSTMap has been developed at Kurtzman Lab Lehman College
#         For more details, please visit
#         sstmap.org @ https://github.com/KurtzmanLab/SSTMap
#         """
#     classifications = [["SSTMap Analysis", "SSTMap HSA"]]
#
#     tags = [tag for lists in classifications for tag in lists]
#
#     # Override defaults for some parameters
#     parameter_overrides = {
#         "memory_mb": {"default": 6000},
#         "spot_policy": {"default": "Allowed"},
#         "prefetch_count": {"default": 1},  # 1 molecule at a time
#         "item_count": {"default": 1}  # 1 molecule at a time
#     }
#
#     start_frame = parameter.IntegerParameter(
#         'start_frame',
#         default=0,
#         min_value=0,
#         level=ADVANCED,
#         help_text="Frame index to start the SSTMap analysis. Default: 0."
#     )
#
#     total_frames = parameter.IntegerParameter(
#         'total_frames',
#         max_value=100000,
#         default=100,
#         level=ADVANCED,
#         help_text="Total number of frames to process during the analysis. Default: 100."
#     )
#
#     hsa_rad = parameter.DecimalParameter(
#         'hsa_rad',
#         min_value=5.0,
#         max_value=10.0,
#         default=5.0,
#         help_text="Distance cutoff (in Angstrom) used to identify hsa region. All waters within this distance from any"
#                   " of the ligand atom are included in the analysis. Default: 5.0."
#     )
#
#     lig_res_name = parameter.StringParameter(
#         'lig_res_name',
#         default='LIG',
#         max_length=4,
#         help_text="Resname to use to identify the ligand"
#     )
#
#     # Values taken from AMBER 17 manual p. 610
#     # Water Model    Mean Energy (Eww-norm) (kcal/mol/water)    Number Density (A^-3)
#     # TIP3P          -9.533                                     0.0329
#     # TIP4PEW        -11.036                                    0.0332
#     # TIP4P          -9.856                                     0.0332
#     # TIP5P          -9.596                                     0.0329
#     # TIP3PFW        -11.369                                    0.0334
#     # SPCE           -11.123                                    0.0333
#     # SPCFW          -11.873                                    0.0329
#     # OPC                                                       0.0333
#
#     wat_model=parameter.StringParameter(
#         'wat_model',
#         default='TIP3P',
#         choices=['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC'],
#         level=ADVANCED,
#         help_text="Water model used during the simulation. Used to set bulk density number. Default: TIP3P."
#     )
#
#     ligand_port = RecordInputPort("ligand_port", initializer=True)
#
#     # Uncomment this and implement if you need to initialize the cube
#     def begin(self):
#         self.opt = vars(self.args)
#         self.opt['Logger'] = self.log
#
#         # Generate dictionary of water models and bulk density
#         wat_model = ['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC']
#         wat_model_bulk_density = [0.0329, 0.0332, 0.0332, 0.0329, 0.0334, 0.0333, 0.0329, 0.0333]
#         self.wat_model_density_dic = dict(zip(wat_model,wat_model_bulk_density))
#
#         for ligand in self.ligand_port:
#             self.ligand = ligand.get_value(Fields.primary_molecule)
#
#     # Records are passed to this function for processing.
#     def process(self, record, port):
#         try:
#             if port == "intake":
#
#                 # Generation options dictionary
#                 opt = dict(self.opt)
#
#                 # Assigning the bulk density bas on water model
#                 opt['rho_bulk'] = self.wat_model_density_dic[opt['wat_model']]
#
#                 if not record.has_value(Fields.primary_molecule):
#                     self.log.error("Missing molecule Primary Molecule' field")
#                     self.failure.emit(record)
#                     return
#
#                 system = record.get_value(Fields.primary_molecule)
#
#                 if not record.has_value(Fields.title):
#                     self.log.warn("Missing record Title field")
#                     system_title = system.GetTitle()[0:12]
#                 else:
#                     system_title = record.get_value(Fields.title)
#
#                 if not record.has_value(Fields.id):
#                     raise ValueError("Missing ID Field")
#
#                 sys_id = record.get_value(Fields.id)
#
#                 sys_info = system_title + '_' + str(sys_id)
#
#                 # Get the MDStageRecord list from the record
#                 if record.has_value(Fields.md_stages):
#                     mdstages = record.get_value(Fields.md_stages)
#                 else:
#                     raise ValueError("Field md_stages is missing!")
#
#                 # Get the MDStageRecord for the production stage from the MDStageRecord list
#                 # That correspond to the last member of the MDStageRecord list
#                 mdstage_prod = mdstages[-1]
#
#                 # Get the MDSystemRecord for the production stage
#                 if mdstage_prod.has_value(Fields.md_system):
#                     mdsystem_prod = mdstage_prod.get_value(Fields.md_system)
#                 else:
#                     raise ValueError("Field md_system is missing!")
#
#                 # Get the PARMED object from the MDSystemRecord
#                 prod_topology_parmed = mdsystem_prod.get_value(Fields.structure)
#
#                 # Generate the parameter supplementary file using the parmed object
#                 with TemporaryDirectory() as output_directory:
#                     opt['Logger'].info("{} - Temporary directory: {}\n".format(sys_info, output_directory))
#
#                     # Get the name of the trajectory from the  production MDStageRecord
#                     if mdstage_prod.has_value(Fields.trajectory):
#
#                         if in_orion():
#                             prod_traj_path = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                             prod_traj_filename = os.path.join(output_directory, "trajectory.h5")
#                             copyfile(prod_traj_path, prod_traj_filename)
#                         else:
#                             prod_traj_filename = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                     else:
#                         raise ValueError("MD_stages do not have a trajectory!")
#
#                     opt['Logger'].info("{} - Trajectory file name: {}".format(sys_info, prod_traj_filename))
#
#                     # Get the final structure from the production stage
#                     prod_coord_eomol = mdsystem_prod.get_value(Fields.topology)
#
#                     # Extract the ligand from the final frame
#                     prot, lig, wat, excp = split(prod_coord_eomol)
#
#                     self.log.info("System name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
#                                   "Water atom numbers = {}\nExcipients atom numbers = {}".format(sys_info,
#                                                                                                  prot.NumAtoms(),
#                                                                                                  lig.NumAtoms(),
#                                                                                                  wat.NumAtoms(),
#                                                                                                  excp.NumAtoms()))
#
#                     # Generate Variables needed for running SSTMap
#                     ligand_filename = os.path.join(output_directory, "ligand.pdb")
#
#                     ofs = oechem.oemolostream(ligand_filename)
#                     pdb_flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
#                     ofs.SetFlavor(oechem.OEFormat_PDB, pdb_flavor)
#                     if lig.GetMaxAtomIdx() > 0:
#                         oechem.OEWriteConstMolecule(ofs, lig)
#                         ligand_align = False
#                     else:
#                         oechem.OEWriteConstMolecule(ofs, self.ligand)
#                         ligand_align = True
#                     ofs.close()
#
#                     if ligand_align:
#                         mdstage_setup = mdstages[0]
#                         if mdstage_setup.get_value(Fields.stage_type) == MDStageTypes.SETUP:
#                             mdsytem_setup = mdstage_setup.get_value(Fields.md_system)
#                             setup_topology = mdsytem_setup.get_value(Fields.topology)
#
#                     opt['Logger'].info("{} - Processing Trajectory".format(sys_info))
#                     # In order to get meaningful energy values we need to strip the ions from the trajectory
#                     # Also, we need to fit the protein to a reference frame to remove translations and
#                     # rotations.
#                     top_filename, parm_filename, prod_traj_filename, aligned_prot_oemol = sstmap_utils.\
#                         process_trajectory(prod_traj_filename,
#                                            prod_topology_parmed,
#                                            opt['lig_res_name'],
#                                            output_directory,
#                                            reference_topology=setup_topology)
#
#                     # Get number of frames
#                     total_number_frames = 0
#                     if prod_traj_filename.endswith(".h5"):
#                         for chunk in md.iterload(prod_traj_filename):
#                             total_number_frames += chunk.n_frames
#                     else:
#                         for chunk in md.iterload(prod_traj_filename, top=top_filename):
#                             total_number_frames += chunk.n_frames
#
#                     # Start SSTMap HSA analysis
#                     # Change to tmp directory to avoid data overwrite
#                     cwd = os.getcwd()
#                     os.chdir(output_directory)
#
#                     opt['Logger'].info("{} - Starting HSA Calculation....".format(sys_info))
#
#                     # Initialize HSA calculation
#                     hsa = sm.SiteWaterAnalysis(topology_file=top_filename,
#                                                trajectory=prod_traj_filename,
#                                                start_frame=0,
#                                                num_frames=total_number_frames,
#                                                supporting_file=parm_filename,
#                                                ligand_file=ligand_filename,
#                                                hsa_region_radius=opt['hsa_rad'])
#
#                     # Initialize hydration sites
#                     hsa.initialize_hydration_sites()
#
#                     # Print System summary
#                     hsa.print_system_summary()
#
#                     # Get frame information
#                     cluster_frame_info_list = cp.deepcopy(hsa.site_waters)
#
#                     # Generate clusters and calculate quantities
#                     hsa.calculate_site_quantities()
#
#                     # Write Calculation summary
#                     hsa.write_calculation_summary()
#
#                     # Write data
#                     hsa.write_data()
#                     ############
#
#                     # Generate EOMol for Water cluster depiction
#                     multi_confomer_cluster_list = sstmap_utils.process_clusters(output_directory, total_number_frames, cluster_frame_info_list)
#
#                     # Generate OEMOL for Most probable configuration
#                     most_prob_config = sstmap_utils.probable_conf(output_directory)
#
#                     # Create new record with results
#                     new_record = OERecord()
#                     new_record.set_value(Fields.primary_molecule, most_prob_config)
#
#                     new_record.set_value(Fields.protein, aligned_prot_oemol)
#
#                     # Create Field
#                     mol_vec = OEField("clusters", Types.Chem.MolVec)
#                     new_record.set_value(mol_vec, multi_confomer_cluster_list)
#
#                     alloutfile = os.path.join(output_directory, 'MC_all_clusters.oeb')
#                     allofs = oechem.oemolostream(alloutfile)
#
#                     for conf in multi_confomer_cluster_list:
#                         oechem.OEWriteConstMolecule(allofs, conf)
#
#                     allofs.close()
#
#                     hsa_data = os.path.join(cwd, "HSA_Results_data")
#                     shutil.copytree(output_directory, hsa_data)
#
#                     self.success.emit(new_record)
#
#         except:
#             # Attach an error message to the molecule that failed
#             self.log.error(traceback.format_exc())
#             # Return failed mol
#             self.failure.emit(record)
#
#
# class SSTMapGist(ParallelMixin, OERecordComputeCube):
#     # Cube documentation.  This documentation for this cube, and all other cubes in this repository, can be converted
#     # to html by calling 'invoke docs' from the root directory of this repository.  This documentation will also
#     # appear in the Orion Floe editor.
#     # version = "0.1.4"
#
#     title = "SSTMAP GIST Analysis"
#
#     description = """
#         SSTMap performs Water Thermodynamics analysis.
#         SSTMaps supports hydration site analysis (HSA)
#         and Grid Inhomogeneous Solvation Theory (GIST).
#
#         SSTMap has been developed at Kurtzman Lab Lehman College
#         For more details, please visit
#         sstmap.org @ https://github.com/KurtzmanLab/SSTMap
#         """
#     classifications = [["SSTMap Analysis", "GIST"]]
#
#     tags = [tag for lists in classifications for tag in lists]
#
#     # Override defaults for some parameters
#     parameter_overrides = {
#         "memory_mb": {"default": 6000},
#         "spot_policy": {"default": "Allowed"},
#         "prefetch_count": {"default": 1},  # 1 molecule at a time
#         "item_count": {"default": 1}  # 1 molecule at a time
#     }
#
#     # The first variable passed to a parameter must always be the variable the parameter is assigned to as a string.
#     grid_res = parameter.DecimalParameter(
#         'grid_res',
#         default=0.5,
#         max_value=0.75,
#         min_value=0.2,
#         level=ADVANCED,
#         help_text='Grid resolution in A. Default: 0.5.'
#     )
#
#     grid_dim = parameter.IntegerParameter(
#         'grid_dim',
#         default=48,
#         level=ADVANCED,
#         help_text="Number of voxels in each direction. Usually grids are square, All dimensions are the same."
#                   "Default: 48."
#     )
#
#     wat_model = parameter.StringParameter(
#         'wat_model',
#         default='TIP3P',
#         choices=['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC'],
#         level=ADVANCED,
#         help_text="Water model used during the simulation. Used to set bulk density number. Default: TIP3P."
#     )
#
#     lig_res_name = parameter.StringParameter(
#         'lig_res_name',
#         default='LIG',
#         max_length=4,
#         help_text="Resname to use to identify the ligand"
#     )
#
#     start_frame = parameter.IntegerParameter(
#         'start_frame',
#         default=0,
#         min_value=0,
#         level=ADVANCED,
#         help_text="Frame index to start the SSTMap analysis. Default: 0."
#
#     )
#
#     total_frames = parameter.IntegerParameter(
#         'total_frames',
#         max_value=100000,
#         default=100,
#         level=ADVANCED,
#         help_text="Total number of frames to process during the analysis. Default: 100."
#     )
#
#     ligand_port = RecordInputPort("ligand_port", initializer=True)
#
#     def begin(self):
#         self.opt = vars(self.args)
#         self.opt['Logger'] = self.log
#
#         # Generate dictionary of water models and bulk density
#         wat_model = ['TIP3P', 'TIP4PEW', 'TIP4P', 'TIP5P', 'TIP3PFW', 'SPCE', 'SPCFW', 'OPC']
#         wat_model_bulk_density = [0.0329, 0.0332, 0.0332, 0.0329, 0.0334, 0.0333, 0.0329, 0.0333]
#         self.wat_model_density_dic = dict(zip(wat_model, wat_model_bulk_density))
#
#         for ligand in self.ligand_port:
#             self.ligand = ligand.get_value(Fields.primary_molecule)
#
#     # Records are passed to this function for processing.
#     def process(self, record, port):
#         try:
#             if port == "intake":
#
#                 # Generation options dictionary
#                 opt = dict(self.opt)
#
#                 # Assigning the bulk density bas on water model
#                 opt['rho_bulk'] = self.wat_model_density_dic[opt['wat_model']]
#
#                 if not record.has_value(Fields.primary_molecule):
#                     self.log.error("Missing molecule Primary Molecule' field")
#                     self.failure.emit(record)
#                     return
#
#                 system = record.get_value(Fields.primary_molecule)
#
#                 if not record.has_value(Fields.title):
#                     self.log.warn("Missing record Title field")
#                     system_title = system.GetTitle()[0:12]
#                 else:
#                     system_title = record.get_value(Fields.title)
#
#                 if not record.has_value(Fields.id):
#                     raise ValueError("Missing ID Field")
#
#                 sys_id = record.get_value(Fields.id)
#
#                 sys_info = system_title + '_' + str(sys_id)
#
#                 # Get the MDStageRecord list from the record
#                 if record.has_value(Fields.md_stages):
#                     mdstages = record.get_value(Fields.md_stages)
#                 else:
#                     raise ValueError("Field md_stages is missing!")
#
#                 # Get the MDStageRecord for the production stage from the MDStageRecord list
#                 # That correspond to the last member of the MDStageRecord list
#                 mdstage_prod = mdstages[-1]
#
#                 # Get the MDSystemRecord for the production stage
#                 if mdstage_prod.has_value(Fields.md_system):
#                     mdsystem_prod = mdstage_prod.get_value(Fields.md_system)
#                 else:
#                     raise ValueError("Field md_system is missing!")
#
#                 # Get the PARMED object from the MDSystemRecord
#                 prod_topology_parmed = mdsystem_prod.get_value(Fields.structure)
#
#                 with TemporaryDirectory() as output_directory:
#                     opt['Logger'].info("{} - Temporary directory: {}\n".format(sys_info , output_directory))
#
#                     # Get the name of the trajectory from the  production MDStageRecord
#                     if mdstage_prod.has_value(Fields.trajectory):
#
#                         if in_orion():
#                             prod_traj_path = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                             prod_traj_filename = os.path.join(output_directory, "trajectory.h5")
#                             copyfile(prod_traj_path, prod_traj_filename)
#                         else:
#                             prod_traj_filename = omm_utils.download(mdstage_prod.get_value(Fields.trajectory))
#                     else:
#                         raise ValueError("MD_stages do not have a trajectory!")
#
#                     opt['Logger'].warn("{} - Trajectory file name: {}".format(sys_info, prod_traj_filename))
#
#                     # Get the final structure from the production stage
#                     prod_coord_eomol = mdsystem_prod.get_value(Fields.topology)
#
#                     # Extract the ligand from the final frame
#                     prot, lig, wat, excp = split(prod_coord_eomol)
#
#                     self.log.info("System name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
#                                   "Water atom numbers = {}\nExcipients atom numbers = {}".format(sys_info,
#                                                                                                  prot.NumAtoms(),
#                                                                                                  lig.NumAtoms(),
#                                                                                                  wat.NumAtoms(),
#                                                                                                  excp.NumAtoms()))
#                     # Generate Variables needed for running SSTMap
#                     ligand_filename = os.path.join(output_directory, "ligand.pdb")
#
#                     opt['Logger'].info("{} - Generate files for GIST....".format(sys_info))
#
#                     # Generate files for GIST
#                     ofs = oechem.oemolostream(ligand_filename)
#
#                     pdb_flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
#
#                     ofs.SetFlavor(oechem.OEFormat_PDB, pdb_flavor)
#
#                     if lig.GetMaxAtomIdx() > 0:
#
#                         oechem.OEWriteConstMolecule(ofs, lig)
#
#                         ligand_align = False
#                     else:
#
#                         oechem.OEWriteConstMolecule(ofs, self.ligand)
#
#                         ligand_align = True
#                     ofs.close()
#
#                     if ligand_align:
#
#                         mdstage_setup = mdstages[0]
#
#                         if mdstage_setup.get_value(Fields.stage_type) == MDStageTypes.SETUP:
#                             mdsytem_setup = mdstage_setup.get_value(Fields.md_system)
#                             setup_topology = mdsytem_setup.get_value(Fields.topology)
#
#                     opt['Logger'].info("{} - Processing Trajectory".format(sys_info))
#
#                     # In order to get meaningful energy values we need to strip the ions from the trajectory
#                     # Also, we need to fit the protein to a reference frame to remove translations and
#                     # rotations.
#                     top_filename, parm_filename, prod_traj_filename, aligned_prot_oemol = sstmap_utils.\
#                         process_trajectory(prod_traj_filename, prod_topology_parmed,
#                                            opt['lig_res_name'], output_directory,
#                                            reference_topology=setup_topology)
#                     # Get number of frames
#                     total_number_frames = 0
#                     if prod_traj_filename.endswith(".h5"):
#                         for chunk in md.iterload(prod_traj_filename):
#                             total_number_frames += chunk.n_frames
#                     else:
#                         for chunk in md.iterload(prod_traj_filename, top=top_filename):
#                             total_number_frames += chunk.n_frames
#
#                     # Start SSTMap GIST analysis
#                     # Change to tmp directory to avoid data overwrite
#                     cwd = os.getcwd()
#                     os.chdir(output_directory)
#
#                     opt['Logger'].info("{} - Starting GIST Calculation".format(sys_info))
#
#                     # Initialize GIST calculation
#                     gist = sm.GridWaterAnalysis(topology_file=top_filename,
#                                                 trajectory=prod_traj_filename,
#                                                 start_frame=0,
#                                                 num_frames=total_number_frames,
#                                                 supporting_file=parm_filename,
#                                                 ligand_file=ligand_filename,
#                                                 grid_dimensions=[opt['grid_dim'], opt['grid_dim'], opt['grid_dim']],
#                                                 grid_resolution=[opt['grid_res'], opt['grid_res'], opt['grid_res']])
#
#                     # Create new record with results
#                     new_record = OERecord()
#                     new_record.set_value(Fields.primary_molecule, aligned_prot_oemol)
#                     new_record.set_value(Fields.title, system_title),
#                     new_record.set_value(Fields.id, sys_id)
#
#                     # Print System summary from GISt
#                     gist.print_system_summary()
#
#                     # Make GIST calculations
#                     gist.calculate_grid_quantities(hbonds=True)
#
#                     # Write GIST Data
#                     gist.write_data()
#
#                     # Generate constat to remove the density weight
#                     g_const = opt['grid_res'] * opt['grid_res'] * opt['grid_res']
#
#                     # Extract voxel coordinates
#                     x_gist_voxels = gist.voxeldata[:, GISTFields.x]
#                     y_gist_voxels = gist.voxeldata[:, GISTFields.y]
#                     z_gist_voxels = gist.voxeldata[:, GISTFields.z]
#
#                     # Extract the data use to generate the OEGrids
#                     g_gO = gist.voxeldata[:, GISTFields.g_O]
#                     g_gH = gist.voxeldata[:, GISTFields.g_H]
#                     g_Eww = gist.voxeldata[:, GISTFields.E_ww_dens]
#                     g_Esw = gist.voxeldata[:, GISTFields.E_sw_dens]
#                     g_So = gist.voxeldata[:, GISTFields.TS_or_dens]
#                     g_St = gist.voxeldata[:, GISTFields.TS_tr_dens]
#
#                     # Remove the density weight
#                     for g_data in [g_St, g_So, g_Eww, g_Esw]:
#                         g_data = g_data * g_const
#
#                     # Calculate total Energy = Ewat-wat + Esolute-wat
#                     g_Etot = g_Esw + g_Eww
#
#                     # Calculate total Entropy = Sorient + Strans
#                     g_Stot = g_So + g_St
#
#                     # Calculate Helmholtz free energy
#                     g_A = g_Etot - g_Stot
#
#                     # create dict for calc values
#                     calc_values = {97: 'Etot',
#                                    98: 'Stot',
#                                    99: 'FreeE'}
#
#                     # Write data to OEGrid
#                     grid_data_field_num = [GISTFields.g_O, GISTFields.g_H,
#                                            GISTFields.E_ww_dens, GISTFields.E_sw_dens,
#                                            GISTFields.TS_or_dens, GISTFields.TS_tr_dens, 97, 98, 99]
#
#                     grid_data = [g_gO, g_gH, g_Eww, g_Esw, g_So, g_St, g_Etot, g_Stot, g_A]
#
#                     grid_data_comb = list(zip(grid_data, grid_data_field_num))
#
#                     # Get grid center
#                     grid_center = gist.center.tolist()
#                     grid_orig = gist.origin.tolist()
#
#                     # Write grid information necessary to recreate the grids from the summary
#                     grid_info_fn = os.path.join(output_directory, "gist_grid_data.txt")
#
#                     with open(grid_info_fn, "w") as g_ofs:
#                         g_ofs.write("grid dimensions: {} {} {}\ngrid center: {} {} {}\ngrid origin: {} {} {}\n"
#                                     "grid resolution: {}\n".format(opt['grid_dim'], opt['grid_dim'], opt['grid_dim'],
#                                                                    grid_center[0], grid_center[1], grid_center[2],
#                                                                    grid_orig[0], grid_orig[1], grid_orig[2],
#                                                                    opt['grid_res']))
#                     for data, data_num in grid_data_comb:
#
#                         if data_num < 90:
#                             grid_name = GISTFields.data_titles[data_num]
#                         else:
#                             grid_name = calc_values[data_num]
#
#                         # Initializing the OEGrid object
#                         grid = oegrid.OEScalarGrid(opt['grid_dim'], opt['grid_dim'], opt['grid_dim'],
#                                                    grid_center[0],
#                                                    grid_center[1], grid_center[2], opt['grid_res'])
#
#                         grid.SetTitle(grid_name)
#
#                         # Set values in OEGrids
#                         for data_pnt in range(grid.GetSize()):
#                             x = x_gist_voxels[data_pnt]
#                             y = y_gist_voxels[data_pnt]
#                             z = z_gist_voxels[data_pnt]
#                             grid.SetValue(x, y, z, data[data_pnt])
#
#                         grid_field = OEField(grid_name, Types.Chem.Grid)
#                         new_record.set_value(grid_field, grid)
#
#                         # Write OEGrid
#                         grid_file_name = grid_name + ".grd"
#                         grid_file_name_wpath = os.path.join(output_directory, grid_file_name)
#                         oegrid.OEWriteGrid(grid_file_name_wpath, grid)
#
#                     # Write dx file of all calculated quantities
#                     gist.generate_dx_files()
#
#                     # Print Calculation Summary
#                     gist.print_calcs_summary()
#
#                     gist_data = os.path.join(cwd, "GIST_Results_data")
#                     shutil.copytree(output_directory, gist_data)
#
#                     self.success.emit(new_record)
#
#         except:
#             # Attach an error message to the molecule that failed
#             self.log.error(traceback.format_exc())
#             print(traceback.format_exc(), flush=True)
#             # Return failed mol
#             self.failure.emit(record)
