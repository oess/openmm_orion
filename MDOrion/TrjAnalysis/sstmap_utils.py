# from openeye import oechem
# import mdtraj as md
# import os
# import glob
#
#
# def process_trajectory(trajectory_filename, parmed_obj, lig_resname, output_dir, reference_topology=None):
#
#     # Save topology
#     if trajectory_filename.endswith(".h5"):
#         topology_filename = os.path.join(output_dir, "strip_topology.h5")
#         top_ifs = open(topology_filename, "w")
#         top_ifs.write("\n")
#         top_ifs.close()
#         parm_filename = os.path.join(output_dir, "parms.txth5")
#         strip_topology_filename = topology_filename
#     else:
#         topology_filename = os.path.join(output_dir, "topology.pdb")
#         parmed_obj.save(topology_filename)
#         # Strip ions from parmed object and save topology
#         strip_topology_filename = os.path.join(output_dir, "strip_topology.pdb")
#         parm_filename = os.path.join(output_dir, "parms.txt")
#
#     parmed_obj.strip('@NA,CL')
#
#     if not trajectory_filename.endswith(".h5"):
#         parmed_obj.save(strip_topology_filename)
#
#     # Generate supporting parameter file
#     with open(parm_filename, "w") as ofs:
#         natoms = len(parmed_obj.atoms)
#         for at in range(natoms):
#             if parm_filename.endswith(".txth5"):
#                 ofs.write("{} {} {}\n".format(parmed_obj.atoms[at].charge, parmed_obj.atoms[at].sigma,
#                                               parmed_obj.atoms[at].epsilon))
#             else:
#                 ofs.write("{} {} {}\n".format(parmed_obj.atoms[at].charge * 18.2223, parmed_obj.atoms[at].sigma,
#                                               parmed_obj.atoms[at].epsilon))
#
#     # strip ions from trajectory
#     if trajectory_filename.endswith(".h5"):
#         trajectory = md.load(trajectory_filename)
#     else:
#         trajectory = md.load(trajectory_filename, top=topology_filename)
#
#     topology_fromtrajectory = trajectory.topology
#     atoms_tokeep = topology_fromtrajectory.select("water or protein or resn ACE or resn NMA or resn NME or resn {}".format(lig_resname))
#     atoms_tofit = topology_fromtrajectory.select("protein and name CA")
#
#     # Fit trajectory to first frame if has ligand or to a reference
#     if reference_topology is None:
#         trajectory_fit = trajectory.superpose(trajectory, frame=0, atom_indices=atoms_tofit, parallel=False)
#         protein_selection = topology_fromtrajectory.select("protein")
#     else:
#         # Generate pdb file for reference topology
#         reference_topology_filename = os.path.join(output_dir, "reference_topology.pdb")
#         ofs = oechem.oemolostream(reference_topology_filename)
#         pdb_flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
#         ofs.SetFlavor(oechem.OEFormat_PDB, pdb_flavor)
#         oechem.OEWriteConstMolecule(ofs, reference_topology)
#         ofs.close()
#
#         # Create reference with mdtraj
#         reference_topology_mdtraj = md.load_pdb(reference_topology_filename)
#         ref_atoms_tofit = reference_topology_mdtraj.topology.select("protein and name CA")
#         protein_selection = reference_topology_mdtraj.topology.select("protein")
#         print("traj: {}, frm: {}".format(len(atoms_tofit), len(ref_atoms_tofit)))
#         trajectory_fit = trajectory.superpose(reference_topology_mdtraj, frame=0,
#                                               atom_indices=ref_atoms_tofit,
#                                               ref_atom_indices=ref_atoms_tofit,
#                                               parallel=False)
#
#     trajectory_fit.atom_slice(atoms_tokeep, inplace=True)
#
#     # Save strip and fit trajectory
#     if trajectory_filename.endswith(".h5"):
#         strip_trajectory_filename = os.path.join(output_dir,"stripandfit_trajectory.h5")
#         trajectory_fit.save_hdf5(filename=strip_trajectory_filename)
#     else:
#         strip_trajectory_filename = os.path.join(output_dir,"stripandfit_trajectory.dcd")
#         trajectory_fit.save_dcd(filename=strip_trajectory_filename)
#
#     # Generate protein for visualization
#     aligned_prot_filename = os.path.join(output_dir, "aligned_protein.pdb")
#     traj_frame = trajectory_fit[1]
#     traj_frame.atom_slice(protein_selection, inplace=True)
#     traj_frame.save_pdb(aligned_prot_filename)
#     # Create protein oemol
#     aligned_prot_oemol = oechem.OEMol()
#     ifs = oechem.oemolistream(aligned_prot_filename)
#     pdb_flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
#     ifs.SetFlavor(oechem.OEFormat_PDB, pdb_flavor)
#     oechem.OEReadMolecule(ifs, aligned_prot_oemol)
#
#     return strip_topology_filename, parm_filename, strip_trajectory_filename, aligned_prot_oemol
#
#
# def process_clusters(output_dir, total_number_frames, cluster_frame_info_list):
#
#     # Change directory to the tmp directory
#     os.chdir(output_dir)
#
#     # Get list of file
#     clusters_list = glob.glob('cluster.*.pdb')
#
#     # Sort the list
#     clusters_list.sort()
#
#     # Dummy coord for water
#     dummy_wat_filename = gen_wat_dummy_file(output_dir)
#     tmp_dummy_wat_mol = oechem.OEMol()
#     ifs_dummy_wat = oechem.oemolistream()
#     ifs_dummy_wat.SetFormat(oechem.OEFormat_PDB)
#     ifs_dummy_wat.open(dummy_wat_filename)
#     oechem.OEReadMolecule(ifs_dummy_wat, tmp_dummy_wat_mol)
#     ifs_dummy_wat.close()
#     hv = oechem.OEHierView(tmp_dummy_wat_mol)
#     for chain in hv.GetChains():
#         for frag in chain.GetFragments():
#             for hres in frag.GetResidues():
#                 hres_atoms = hres.GetAtoms()
#                 atom_members = oechem.OEIsAtomMember(hres_atoms)
#                 dummy_wat_mol = oechem.OEMol()
#                 oechem.OESubsetMol(dummy_wat_mol, tmp_dummy_wat_mol, atom_members)
#
#     # Create list to storage the cluster confomers
#     multi_confomer_cluster_list = []
#
#     # Loop over all the clusters to create
#     # the multi conformer oeb file
#     for clt_num, clt_file in enumerate(clusters_list):
#         # Setup file stream
#         ifs_clts = oechem.oemolistream()
#         ifs_clts.SetFormat(oechem.OEFormat_PDB)
#         ifs_clts.open(clt_file)
#
#         # Read the file
#         mol_clt = oechem.OEMol()
#         oechem.OEReadMolecule(ifs_clts, mol_clt)
#
#         # Create a hierarchical view for easy selection of residues
#         hv = oechem.OEHierView(mol_clt)
#
#         # Loop to select residue
#         ref_frame = 0
#         for chain in hv.GetChains():
#             for frag in chain.GetFragments():
#                 for hres in frag.GetResidues():
#
#                     # Make a selection for the residue
#                     hres_atoms = hres.GetAtoms()
#                     atom_members = oechem.OEIsAtomMember(hres_atoms)
#                     conf = oechem.OEMol()
#                     oechem.OESubsetMol(conf, mol_clt, atom_members)
#
#                     # Initialize the container with the firs conformer
#                     if hres.GetResidueNumber() == 0 and cluster_frame_info_list[clt_num][hres.GetResidueNumber()][0]\
#                             == 0 and ref_frame == 0:
#                         multy_conf_cluster = oechem.OEMol(conf)
#                         ref_frame += 1
#                         continue
#                     elif hres.GetResidueNumber() == 0 and cluster_frame_info_list[clt_num][hres.GetResidueNumber()][0]\
#                             != 0 and ref_frame == 0:
#                         multy_conf_cluster = oechem.OEMol(dummy_wat_mol)
#                         ref_frame += 1
#                         for i in range(1, cluster_frame_info_list[clt_num][hres.GetResidueNumber()][0]):
#                             multy_conf_cluster.NewConf(dummy_wat_mol)
#                             ref_frame += 1
#
#                         multy_conf_cluster.NewConf(conf)
#                         ref_frame += 1
#                         continue
#
#                     # If we are in the correct frame then add conformer if not add dummy
#                     if ref_frame == cluster_frame_info_list[clt_num][hres.GetResidueNumber()][0]:
#                         # Add conformer
#                         multy_conf_cluster.NewConf(conf)
#                         ref_frame += 1
#                         continue
#                     else:
#                         for i in range(ref_frame, cluster_frame_info_list[clt_num][hres.GetResidueNumber()][0]):
#                             multy_conf_cluster.NewConf(dummy_wat_mol)
#                             ref_frame += 1
#
#                         multy_conf_cluster.NewConf(conf)
#                         ref_frame += 1
#                         continue
#
#         if ref_frame < total_number_frames:
#             for n in range(ref_frame, total_number_frames):
#                 multy_conf_cluster.NewConf(dummy_wat_mol)
#                 ref_frame += 1
#
#         # Print OEMol
#         num_conf = multy_conf_cluster.GetMaxConfIdx()
#         if not num_conf < total_number_frames:
#             # Save all multi conformers in a single list and return this object
#             multi_confomer_cluster_list.append(multy_conf_cluster)
#         else:
#             print(
#                 "Number of conformers ({}) not equal to Number of frames ({})".format(num_conf, total_number_frames))
#
#     return multi_confomer_cluster_list
#
#
# def probable_conf(output_dir):
#
#     # Change to tmp directory
#     os.chdir(output_dir)
#     prob_config_filename = os.path.join(output_dir, "probable_configs.pdb")
#     ifs_prob_config = oechem.oemolistream()
#     ifs_prob_config.SetFormat(oechem.OEFormat_PDB)
#     ifs_prob_config.open(prob_config_filename)
#     prob_config_oemol = oechem.OEMol()
#     oechem.OEReadMolecule(ifs_prob_config, prob_config_oemol)
#
#     # Create a hierarchical view for easy selection of residues
#     hv = oechem.OEHierView(prob_config_oemol)
#
#     # Loop to select residue
#     multy_conf_probconfig = oechem.OEMol()
#     for chain in hv.GetChains():
#         for frag in chain.GetFragments():
#             for hres in frag.GetResidues():
#
#                 # Make a selection for the residue
#                 hres_atoms = hres.GetAtoms()
#                 atom_members = oechem.OEIsAtomMember(hres_atoms)
#                 conf = oechem.OEMol()
#                 oechem.OESubsetMol(conf, prob_config_oemol, atom_members)
#                 oechem.OEAddMols(multy_conf_probconfig, conf)
#
#     return multy_conf_probconfig
#
#
# def gen_wat_dummy_file(output_dir):
#
#     # Change to tmp directory
#     os.chdir(output_dir)
#     dummy_wat_filename = os.path.join(output_dir, "dummy_water.pdb")
#
#     ofs = open(dummy_wat_filename, "w")
#     dummy_wat_atoms="ATOM      0  O   WAT A   0       0.000   0.000   0.000  0.00  0.00           O\n" \
#                     "ATOM      1  H1  WAT A   0       0.958   0.000   0.000  0.00  0.00           H\n" \
#                     "ATOM      2  H2  WAT A   0      -0.239   0.928   0.000  0.00  0.00           H"
#
#     ofs.write(dummy_wat_atoms)
#     ofs.close()
#
#     return dummy_wat_filename
#
#
# class GISTFields:
#     data_titles = ['index', 'x', 'y', 'z',
#                    'N_wat', 'g_O', 'g_H',
#                    'TS_tr_dens', 'TS_tr_norm',
#                    'TS_or_dens', 'TS_or_norm',
#                    'dTSsix-dens', 'dTSsix_norm',
#                    'E_sw_dens', 'E_sw_norm', 'E_ww_dens', 'Eww_norm',
#                    'E_ww_nbr_dens', 'E_ww_nbr_norm',
#                    'N_nbr_dens', 'N_nbr_norm',
#                    'f_hb_dens', 'f_hb_norm',
#                    'N_hb_sw_dens', 'N_hb_sw_norm', 'N_hb_ww_dens', 'N_hb_ww_norm',
#                    'N_don_sw_dens', 'N_don_sw_norm', 'N_acc_sw_dens', 'N_acc_sw_norm',
#                    'N_don_ww_dens', 'N_don_ww_norm', 'N_acc_ww_dens', 'N_acc_ww_norm']
#
#     # Voxel x y z nwat gO gH dTStr-dens dTStr-norm dTSor-dens dTSor-norm
#     # dTSsix-dens dTSsix-norm Esw-dens Esw-norm Eww-dens Eww-norm Eww-nbr-dens Eww-nbr-norm Nnbr-dens Nnbr-norm
#     # fHB-dens fHB-norm Nhbsw_dens Nhbsw_norm Nhbww_dens Nhbww_norm Ndonsw_dens Ndonsw_norm Naccsw_dens Naccsw_norm
#     # Ndonww_dens Ndonww_norm Naccww_dens Naccww_norm
#     index = 0
#     x = 1
#     y = 2
#     z = 3
#     N_wat = 4
#     g_O = 5
#     g_H = 6
#     TS_tr_dens = 7
#     TS_tr_norm = 8
#     TS_or_dens = 9
#     TS_or_norm = 10
#     dTSsix_dens = 11
#     dTSsix_norm = 12
#     E_sw_dens = 13
#     E_sw_norm = 14
#     E_ww_dens = 15
#     Eww_norm = 16
#     E_ww_nbr_dens = 17
#     E_ww_nbr_norm = 18
#     N_nbr_dens = 19
#     N_nbr_norm = 20
#     f_hb_dens = 21
#     f_hb_norm = 22
#     N_hb_sw_dens = 23
#     N_hb_sw_norm = 24
#     N_hb_ww_dens = 25
#     N_hb_ww_norm = 26
#     N_don_sw_dens = 27
#     N_don_sw_norm = 28
#     N_acc_sw_dens = 29
#     N_acc_sw_norm = 30
#     N_don_ww_dens = 31
#     N_don_ww_norm = 32
#     N_acc_ww_dens = 33
#     N_acc_ww_norm = 34
#
#
# class HSAFields:
#     data_titles = ["index", "x", "y", "z",
#                    "nwat", "occupancy",
#                    "Esw", "EswLJ", "EswElec",
#                    "Eww", "EwwLJ", "EwwElec", "Etot", "Ewwnbr",
#                    "TSsw_trans", "TSsw_orient", "TStot",
#                    "Nnbrs", "Nhbww", "Nhbsw", "Nhbtot",
#                    "f_hb_ww", "f_enc",
#                    "Acc_ww", "Don_ww", "Acc_sw", "Don_sw",
#                    "solute_acceptors", "solute_donors"]
#     index = 0
#     x = 1
#     y = 2
#     z = 3
#     nwat = 4
#     occupancy = 5
#     Esw = 6
#     EswLJ = 7
#     EswElec = 8
#     Eww = 9
#     EwwLJ = 10
#     EwwElec = 11
#     Etot = 12
#     Ewwnbr = 13
#     TSsw_trans = 14
#     TSsw_orient = 15
#     TStot = 16
#     Nnbrs = 17
#     Nhbww = 18
#     Nhbsw = 19
#     Nhbtot = 20
#     f_hb_ww = 21
#     f_enc = 22
#     Acc_ww = 23
#     Don_ww = 24
#     Acc_sw = 25
#     Don_sw = 26
#     solute_acceptors = 27
#     solute_donors = 28