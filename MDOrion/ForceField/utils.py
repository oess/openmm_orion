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

from openeye import oechem

from MDOrion.ForceField.ff_library import ff_library

from MDOrion.ForceField.ffutils import (ParamMolStructure,
                                        parametrize_component,
                                        parametrize_unknown_component,
                                        clean_tags)
import parmed

from oeommtools import utils as oeommutils

from simtk.openmm import app

from openeye import oespruce



class ParametrizeDU:
    def __init__(self, design_unit, opt):

        # Design Unit
        self.du = design_unit
        self.opt = opt

        # Force Field to use: User Defined
        self.other_ff = ff_library.otherff[opt['other_forcefield']]
        self.ligand_ff = ff_library.ligandff[opt['ligand_forcefield']]
        self.protein_ff = ff_library.proteinff[opt['protein_forcefield']]

        # Extended Force Field
        self.protein_extended_ff = ff_library.protein_extended_ff[opt['protein_forcefield']]

        # Force Field to use: Default
        self.counter_ions_ff = ff_library.counter_ionsff['Counter_ions']
        self.metals_ff = ff_library.metals_ff['Metals']
        self.excipients_ff = ff_library.excipients_ff['Excipients']
        self.solvent_ff = ff_library.solventff['Tip3p']
        self.cofactors_ff = ff_library.cofactors_ff['Cofactors']
        self.lipids_ff = ff_library.lipids_ff['Lipids']
        self.nucleics_ff = ff_library.nucleics_ff['Nucleics']

        # What to parametrize
        self.protein = None
        self.ligand = None
        self.other_ligands = None
        self.counter_ions = None
        self.metals = None
        self.excipients = None
        self.solvent = None
        self.cofactors = None
        self.other_cofactors = None
        self.lipids = None
        self.nucleics = None
        self.other_nucleics = None

        # Components found in the current DU
        self.components = dict()

        for pair in self.du.GetTaggedComponents():

            comp_name = pair[0]
            comp_id = self.du.GetComponentID(comp_name)
            comp = pair[1]

            self.components[comp_name] = comp

            if comp_id == oechem.OEDesignUnitComponents_Protein:
                self.protein = comp
            elif comp_id == oechem.OEDesignUnitComponents_Ligand:
                self.ligand = comp
            elif comp_id == oechem.OEDesignUnitComponents_OtherLigands:
                self.other_ligand = comp
            elif comp_id == oechem.OEDesignUnitComponents_CounterIons:
                self.counter_ions = comp
            elif comp_id == oechem.OEDesignUnitComponents_Metals:
                self.metals = comp
            elif comp_id == oechem.OEDesignUnitComponents_Excipients:
                self.excipients = comp
            elif comp_id == oechem.OEDesignUnitComponents_Solvent:
                self.solvent = comp
            elif comp_id == oechem.OEDesignUnitComponents_Cofactors:
                self.cofactors = comp
            elif comp_id == oechem.OEDesignUnitComponents_OtherCofactors:
                self.other_cofactors = comp
            elif comp_id == oechem.OEDesignUnitComponents_Lipids:
                self.lipids = comp
            elif comp_id == oechem.OEDesignUnitComponents_Nucleic:
                self.nucleics = comp
            elif comp_id == oechem.OEDesignUnitComponents_OtherNucleics:
                self.other_nucleics = comp
            else:
                self.opt['Logger'].warn("The following component parametrization is not currently supported: {}".format(comp_name))
        if not self.components:
            raise ValueError("None of the DU components cannot be parametrized")

    @property
    def parametrize_protein(self):
        if self.protein:

            # oechem.OEPerceiveResidues(self.protein, oechem.OEPreserveResInfo_All)

            # for oe_res in oechem.OEGetResidues(self.protein):
            #     print(oe_res.GetName(), oe_res.GetResidueNumber())
            #     for oe_at in oechem.OEGetResidueAtoms(self.protein, oe_res, oechem.OEAssumption_BondedResidue +
            #                         oechem.OEAssumption_ResPerceived +
            #                        oechem.OEAssumption_PDBOrder):
            #         print("\t {}".format(oe_at))

            # OpenMM topology and positions from OEMol
            topology, positions = oeommutils.oemol_to_openmmTop(self.protein)

            # Try to apply the selected FF on the Protein
            forcefield = app.ForceField(self.protein_ff)

            if not forcefield.getUnmatchedResidues(topology):
                omm_protein = forcefield.createSystem(topology, rigidWater=False, constraints=None)
                protein_pmd = parmed.openmm.load_topology(topology, omm_protein, xyz=positions)
                return protein_pmd
            else:
                # Try to apply the extended FF to the Protein
                forcefield = app.ForceField(self.protein_extended_ff)
                unmatched_res_list = forcefield.getUnmatchedResidues(topology)
                if unmatched_res_list:
                    raise ValueError("The following residues cannot be parametrized by the selected FF {}\n{}".
                                     format(self.opt['protein_forcefield'], unmatched_res_list))
                else:
                    omm_protein = forcefield.createSystem(topology, rigidWater=False, constraints=None)
                    protein_pmd = parmed.openmm.load_topology(topology, omm_protein, xyz=positions)
                    return protein_pmd
        else:
            raise ValueError("Protein is not present in the DU")
        pass

    @property
    def parametrize_ligand(self):

        if self.ligand:
            prefix_name = 'LIG'
            pmd = ParamMolStructure(self.ligand, self.ligand_ff, prefix_name=prefix_name)
            ligand_pmd = pmd.parameterize()
            ligand_pmd.residues[0].name = prefix_name
            return ligand_pmd
        else:
            raise ValueError("Ligand is not present in the DU")

    @property
    def parametrize_other_ligands(self):

        if self.other_ligands:
            ligand_pmd = parametrize_unknown_component(self.other_ligands, self.ligand_ff)
            return ligand_pmd
        else:
            raise ValueError("Other Ligands are not present in the DU")

    @property
    def parametrize_counter_ions(self):

        if self.counter_ions:

            counter_ions_pmd, unrec_ions = parametrize_component(self.counter_ions, self.counter_ions_ff)

            if counter_ions_pmd is not None:

                if unrec_ions is not None:

                    unk_ions_pmd = parametrize_unknown_component(unrec_ions, self.other_ff)
                    counter_ions_pmd += unk_ions_pmd

                    self.opt['Logger'].warn("Some Counter Ions have  been parametrized by using the ff: {}".
                                            format(self.other_ff))

                return counter_ions_pmd
            else:
                raise ValueError("Counter Ions cannot be parametrized by using the FF: {}".format(self.counter_ions_ff))
        else:
            raise ValueError("Counter Ions are not present in the DU")

    @property
    def parametrize_metals(self):

        if self.metals:
            metals_pmd, unrec_metals = parametrize_component(self.metals, self.metals_ff)

            if metals_pmd is not None:

                if unrec_metals is not None:

                    unk_metals_pmd = parametrize_unknown_component(unrec_metals, self.other_ff)
                    metals_pmd += unk_metals_pmd

                    self.opt['Logger'].warn("Some Metal Ions have been parametrized by using the ff: {}".
                                            format(self.other_ff))
                return metals_pmd
            else:
                raise ValueError("Metals cannot be parametrized by using the FF: {}".format(self.metals_ff))
        else:
            raise ValueError("Metals are not present in the DU")

    @property
    def parametrize_excipients(self):

        if self.excipients:
            excipients_pmd, unrec_exc = parametrize_component(self.excipients, self.excipients_ff)

            if excipients_pmd is not None:

                if unrec_exc is not None:
                    unk_exc_pmd = parametrize_unknown_component(unrec_exc, self.other_ff)
                    excipients_pmd += unk_exc_pmd

                    self.opt['Logger'].warn("Some Excipients have been parametrized by using the ff: {}".
                                            format(self.other_ff))
                return excipients_pmd
            else:
                raise ValueError("Excipients cannot be parametrized by using the FF: {}".format(self.excipients_ff))
        else:
            raise ValueError("Excipients are not present in the DU")

    @property
    def parametrize_solvent(self):

        if self.solvent:

            oechem.OEPerceiveResidues(self.solvent, oechem.OEPreserveResInfo_All)

            # PDB_Order Error. Call Perceive Residue for the wrong PDB Atom Order
            # for oe_res in oechem.OEGetResidues(self.solvent):
            #     print(oe_res.GetName(), oe_res.GetResidueNumber())
            #     for oe_at in oechem.OEGetResidueAtoms(self.solvent, oe_res, oechem.OEAssumption_BondedResidue +
            #                         oechem.OEAssumption_ResPerceived +
            #                        oechem.OEAssumption_PDBOrder):
            #         print("\t {}".format(oe_at))

            oechem.OEPerceiveResidues(self.solvent, oechem.OEPreserveResInfo_All)

            solvent_pmd, unrec_solvent = parametrize_component(self.solvent, self.solvent_ff)

            if solvent_pmd is not None:

                if unrec_solvent is not None:

                    unk_solvent_pmd = parametrize_unknown_component(unrec_solvent, self.other_ff)
                    solvent_pmd += unk_solvent_pmd

                    self.opt['Logger'].warn("Some Solvent molecules have been parametrized by using the ff: {}".
                                            format(self.other_ff))
                return solvent_pmd
            else:
                raise ValueError("Solvent cannot be parametrized by using the FF: {}".format(self.solvent_ff))
        else:
            raise ValueError("Solvent is not present in the DU")

    @property
    def parametrize_cofactors(self):

        if self.cofactors:
            cofactors_pmd = parmed.Structure()
            unrecognized_cofactors = oechem.OEMol(self.cofactors)

            for idx in range(0, len(self.cofactors_ff)):
                cof_ff = self.cofactors_ff[idx]

                rec_cofactors_pmd, unrecognized_cofactors = parametrize_component(unrecognized_cofactors, cof_ff)

                if rec_cofactors_pmd is not None:
                    cofactors_pmd += rec_cofactors_pmd

            if unrecognized_cofactors is not None:
                unk_cof_pmd = parametrize_unknown_component(unrecognized_cofactors, self.other_ff)
                cofactors_pmd += unk_cof_pmd

            return cofactors_pmd
        else:
            raise ValueError("Cofactors are not present in the DU")

    @property
    def parametrize_other_cofactors(self):

        if self.other_cofactors:
            other_cofactors_pmd = parmed.Structure()
            unrecognized_other_cofactors = oechem.OEMol(self.other_cofactors)

            for idx in range(0, len(self.cofactors_ff)):
                cof_ff = self.cofactors_ff[idx]
                rec_other_cofactors_pmd, unrecognized_other_cofactors = parametrize_component(unrecognized_other_cofactors, cof_ff)

                if rec_other_cofactors_pmd is not None:
                    other_cofactors_pmd += rec_other_cofactors_pmd

                if unrecognized_other_cofactors is None:
                    break

            if unrecognized_other_cofactors is not None:
                unk_other_cof_pmd = parametrize_unknown_component(unrecognized_other_cofactors, self.other_ff)
                other_cofactors_pmd += unk_other_cof_pmd

            return other_cofactors_pmd
        else:
            raise ValueError("Other Cofactors are not present in the DU")

    @property
    def parametrize_lipids(self):

        if self.lipids:
            lipids_pmd, unrec_lipids = parametrize_component(self.lipids, self.lipids_ff)

            if unrec_lipids is not None:
                raise ValueError("Some lipid molecule cannot be parametrized by the FF: {}".format(self.lipids_ff))
            else:
                return lipids_pmd
        else:
            raise ValueError("Lipids are not present in the DU")

    @property
    def parametrize_nucleics(self):

        if self.nucleics:
            nucleics_pmd, unrec_nucleics = parametrize_component(self.nucleics, self.nucleics_ff)

            if unrec_nucleics is not None:
                raise ValueError("Some Nucleics molecule cannot be parametrized by the FF: {}".format(self.nucleics_ff))
            else:
                return nucleics_pmd
        else:
            raise ValueError("Nucleics molecule are not present in the DU")

    @property
    def parametrize_other_nucleics(self):

        if self.nucleics:
            other_nucleics_pmd, other_unrec_nucleics = parametrize_component(self.other_nucleics, self.nucleics_ff)

            if other_unrec_nucleics is not None:
                raise ValueError("Some Other Nucleics molecules cannot be parametrized by the FF: {}".format(self.nucleics_ff))
            else:
                return other_nucleics_pmd
        else:
            raise ValueError("Other Nucleics molecules are not present in the DU")

    @property
    def parametrize_du(self):

        du_pmd = parmed.Structure()
        for comp_name in self.components:

            comp_id = self.du.GetComponentID(comp_name)

            if comp_id == oechem.OEDesignUnitComponents_Protein:
                du_pmd += self.parametrize_protein
            elif comp_id == oechem.OEDesignUnitComponents_Ligand:
                du_pmd += self.parametrize_ligand
            elif comp_id == oechem.OEDesignUnitComponents_OtherLigands:
                du_pmd += self.parametrize_other_ligands
            elif comp_id == oechem.OEDesignUnitComponents_CounterIons:
                du_pmd += self.parametrize_counter_ions
            elif comp_id == oechem.OEDesignUnitComponents_Metals:
                du_pmd += self.parametrize_metals
            elif comp_id == oechem.OEDesignUnitComponents_Excipients:
                du_pmd += self.parametrize_excipients
            elif comp_id == oechem.OEDesignUnitComponents_Solvent:
                du_pmd += self.parametrize_solvent
            elif comp_id == oechem.OEDesignUnitComponents_Cofactors:
                du_pmd += self.parametrize_cofactors
            elif comp_id == oechem.OEDesignUnitComponents_OtherCofactors:
                du_pmd += self.parametrize_other_cofactors
            elif comp_id == oechem.OEDesignUnitComponents_Lipids:
                du_pmd += self.parametrize_lipids
            elif comp_id == oechem.OEDesignUnitComponents_Nucleic:
                du_pmd += self.parametrize_nucleics
            elif comp_id == oechem.OEDesignUnitComponents_OtherNucleics:
                du_pmd += self.parametrize_other_nucleics
            else:
                continue

        return du_pmd


class MDComponents:

    def __init__(self, system_representation, components_title="MD Components"):

        du = None
        molecules = None

        if type(system_representation) == oechem.OEDesignUnit:
            du = system_representation
        elif type(system_representation) == oechem.OEMol:
            molecules = system_representation
        else:
            raise ValueError("The MDComponent class can be initialized with an OE Design Unit "
                             "or an OE Mol. The object passed is: {}".format(type(system_representation)))

        # What to parametrize
        self._protein = None
        self._ligand = None
        self._other_ligands = None
        self._counter_ions = None
        self._metals = None
        self._excipients = None
        self._solvent = None
        self._cofactors = None
        self._other_cofactors = None
        self._lipids = None
        self._nucleics = None
        self._other_nucleics = None

        self._component_names = ['protein', 'ligand',
                                 'other_ligands', 'counter_ions',
                                 'metals', 'excipients',
                                 'solvent', 'cofactors',
                                 'other_cofactors', 'lipids',
                                 'nucleics', 'other_nucleics']

        # Components found in the system_representation
        self._components = dict()

        self._components_title = components_title

        if du is not None:
            self._initialize_from_du(du)
        else:
            self._initialize_from_molecules(molecules)

    def __repr__(self):
        ret_str = "{:<20} {:>7}\n".format("Comp_name", "Atoms")
        ret_str += 28 * "-" + "\n"
        for comp_name, comp in self._components.items():
            ret_str += "{:<20} {:>7}\n".format(comp_name, comp.NumAtoms())
        ret_str += 28 * "-" + "\n"
        ret_str += "{:<20} {:>7}\n".format("Total_Atoms", self.num_atoms)

        return ret_str

    def _initialize_from_du(self, du):

        for pair in du.GetTaggedComponents():

            comp_name = pair[0]
            comp_id = du.GetComponentID(comp_name)
            comp = pair[1]

            # Removing Interaction Hint Container and Style from the components
            oechem.OEDeleteInteractionsHintSerializationData(comp)
            oechem.OEDeleteInteractionsHintSerializationIds(comp)
            oechem.OEClearStyle(comp)

            if comp_id == oechem.OEDesignUnitComponents_Protein:
                self._protein = comp
                self._components['protein'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Ligand:
                self._ligand = comp
                self._components['ligand'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_OtherLigands:
                self._other_ligand = comp
                self._components['other_ligands'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_CounterIons:
                self._counter_ions = comp
                self._components['counter_ions'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Metals:
                self._metals = comp
                self._components['metals'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Excipients:
                self.excipients = comp
                self._components['excipients'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Solvent:
                self._solvent = comp
                self._components['solvent'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Cofactors:
                self._cofactors = comp
                self._components['cofactors'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_OtherCofactors:
                self._other_cofactors = comp
                self._components['other_cofactors'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Lipids:
                self._lipids = comp
                self._components['lipids'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Nucleic:
                self._nucleics = comp
                self._components['nucleics'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_OtherNucleics:
                self._other_nucleics = comp
                self._components['other_nucleics'] = comp
            else:
                print("WARNING: The following component is not currently supported: {}".format(comp_name))
        if not self._components:
            raise ValueError("None of the DU components cannot recognized")

    def _initialize_from_molecules(self, molecules):

        build_opts = oespruce.OEDesignUnitBuildOptions()
        build_opts.SetBuildSidechains(False)
        build_opts.SetBuildLoops(False)
        build_opts.SetCapCTermini(False)
        build_opts.SetCapNTermini(False)
        build_opts.SetDeleteClashingSolvent(False)

        enum_opts = oespruce.OEDesignUnitEnumerateSitesOptions()
        enum_opts.SetAddInteractionHints(False)
        enum_opts.SetAddStyle(False)
        enum_opts.SetEnumerateCofactorSites(False)
        enum_opts.SetDuplicateRemoval(False)
        enum_opts.SetCollapseNonSiteAlts(False)

        prep_opts = oespruce.OEDesignUnitPrepOptions()
        prep_opts.SetProtonate(False)
        prep_opts.SetAssignPartialChargesAndRadii(False)
        prep_opts.SetBuildOptions(build_opts)
        prep_opts.SetEnumerateSitesOptions(enum_opts)

        split_opts = oespruce.OEDesignUnitSplitOptions()
        split_opts.SetMakePackingResidues(False)

        bio_opts = oespruce.OEBioUnitExtractionOptions()
        bio_opts.SetSuperpose(False)

        du_opts = oespruce.OEMakeDesignUnitOptions(split_opts, prep_opts, bio_opts)

        du_meta_data = oespruce.OEStructureMetadata()

        du_list = []

        for du in oespruce.OEMakeDesignUnits(molecules, du_meta_data, du_opts):
            du_list.append(du)

        # Take the first du from the list if available
        if du_list:
            du = du_list[0]
            self._initialize_from_du(du)
        else:
            # Split the complex in components
            protein, ligand, water, excipients = oeommutils.split(molecules,
                                                                  ligand_res_name='LIG')
            solvent = oechem.OEMol()

            if protein.NumAtoms():
                protein = clean_tags(protein)
                self._protein = protein
                self._components['protein'] = protein
            if ligand.NumAtoms():
                ligand = clean_tags(ligand)
                self._ligand = ligand
                self._components['ligand'] = ligand
            if water.NumAtoms():
                water = clean_tags(water)
                oechem.OEAddMols(solvent, water)
            if excipients.NumAtoms():
                excipients = clean_tags(excipients)
                oechem.OEAddMols(solvent, excipients)
            if solvent.NumAtoms():
                self._solvent = solvent
                self._components['solvent'] = solvent

        if self.num_atoms != molecules.NumAtoms():
            raise ValueError("Atom number mismatch: {} vs {}".format(self.num_atoms, molecules.NumAstoms()))

    def __getstate__(self):

        def mol_to_bytes(mol):
            return oechem.OEWriteMolToBytes(oechem.OEFormat_OEB, True, mol)

        state = dict(protein=mol_to_bytes(self._protein) if self._protein else None,
                     ligand=mol_to_bytes(self._ligand) if self._ligand else None,
                     other_ligands=self._other_ligands if self._other_ligands else None,
                     counter_ions=mol_to_bytes(self._counter_ions) if self._counter_ions else None,
                     metals=mol_to_bytes(self._metals) if self._metals else None,
                     excipients=mol_to_bytes(self._excipients) if self._excipients else None,
                     solvent=mol_to_bytes(self._solvent) if self._solvent else None,
                     cofactors=mol_to_bytes(self._cofactors) if self._cofactors else None,
                     other_cofactors=mol_to_bytes(self._other_cofactors) if self._other_cofactors else None,
                     lipids=mol_to_bytes(self._lipids) if self._lipids else None,
                     nucleics=mol_to_bytes(self._nucleics) if self._nucleics else None,
                     other_nucleics=mol_to_bytes(self._other_nucleics) if self._other_nucleics else None,
                     components_title=self._components_title
                     )

        return state

    def __setstate__(self, state):

        def mol_from_bytes(mol_bytes):
            mol = oechem.OEMol()
            oechem.OEReadMolFromBytes(mol, oechem.OEFormat_OEB, True, mol_bytes)
            return mol

        self._components = dict()
        self._total_atoms = 0

        for comp_name, comp in state.items():

            if comp_name == 'components_title':
                self._components_title = comp
                continue

            mol = mol_from_bytes(comp) if comp else None

            if comp_name == 'protein':
                self._protein = mol
            elif comp_name == 'ligand':
                self._ligand = mol
            elif comp_name == 'other_ligands':
                self._other_ligands = mol
            elif comp_name == 'counter_ions':
                self._counter_ions = mol
            elif comp_name == 'metals':
                self._metals = mol
            elif comp_name == 'excipients':
                self._excipients = mol
            elif comp_name == 'solvent':
                self._solvent = mol
            elif comp_name == 'cofactors':
                self._cofactors = mol
            elif comp_name == 'other_cofactors':
                self._other_cofactors = mol
            elif comp_name == 'lipids':
                self._lipids = mol
            elif comp_name == 'nucleics':
                self._nucleics = mol
            elif comp_name == 'other_nucleics':
                self._other_nucleics = mol
            else:
                raise ValueError("Cannot Deserialize Component {}".format(comp_name))

            if mol:
                self._components[comp_name] = mol

    @property
    def num_atoms(self):

        tot_atoms = 0
        for comp_name, comp in self._components.items():
            tot_atoms += comp.NumAtoms()

        return tot_atoms

    @property
    def create_flask(self):
        flask = oechem.OEMol()
        for comp_name, comp in self._components.items():
            if not oechem.OEAddMols(flask, comp):
                raise ValueError("The flask cannot be created. Problems with the component: {}".format(comp_name))
        flask.SetTitle(self._components_title)
        return flask

    @property
    def get_protein(self):
        if self._protein is not None:
            return self._protein
        else:
            raise ValueError("Protein Component has not been found")

    def set_protein(self, protein):
        self._protein = protein

    @property
    def has_protein(self):
        if self._protein is not None:
            return True
        else:
            return False

    @property
    def get_ligand(self):
        if self._ligand is not None:
            return self._ligand
        else:
            raise ValueError("Ligand Component has not been found")

    def set_ligand(self, ligand):
        self._ligand = ligand

    @property
    def has_ligand(self):
        if self._ligand is not None:
            return True
        else:
            return False

    @property
    def get_other_ligands(self):
        if self._other_ligands is not None:
            return self._other_ligands
        else:
            raise ValueError("Other Ligand Component has not been found")

    def set_other_ligands(self, other_ligands):
        self._other_ligands = other_ligands

    @property
    def has_other_ligands(self):
        if self._other_ligands is not None:
            return True
        else:
            return False

    @property
    def get_counter_ions(self):
        if self._counter_ions is not None:
            return self._counter_ions
        else:
            raise ValueError("Counter Ions Component has not been found")

    def set_counter_ions(self, counter_ions):
        self._counter_ions = counter_ions

    @property
    def has_counter_ions(self):
        if self._counter_ions is not None:
            return True
        else:
            return False

    @property
    def get_metals(self):
        if self._metals is not None:
            return self._metals
        else:
            raise ValueError("Metals Component has not been found")

    def set_metals(self, metals):
        self._metals = metals

    @property
    def has_metals(self):
        if self._metals is not None:
            return True
        else:
            return False

    @property
    def get_excipients(self):
        if self._excipients is not None:
            return self._excipients
        else:
            raise ValueError("Excipients Component has not been found")

    def set_excipients(self, excipients):
        self._excipients = excipients

    @property
    def has_excipients(self):
        if self._excipients is not None:
            return True
        else:
            return False

    @property
    def get_solvent(self):
        if self._solvent is not None:
            return self._solvent
        else:
            raise ValueError("Solvent Component has not been found")

    def set_solvent(self, solvent):
        self._solvent = solvent

    @property
    def has_solvent(self):
        if self._solvent is not None:
            return True
        else:
            return False

    @property
    def get_cofactors(self):
        if self._cofactors is not None:
            return self._cofactors
        else:
            raise ValueError("Cofactors Component has not been found")

    def set_cofactors(self, cofactors):
        self._cofactors = cofactors

    @property
    def has_cofactors(self):
        if self._cofactors is not None:
            return True
        else:
            return False

    @property
    def get_other_cofactors(self):
        if self._other_cofactors is not None:
            return self._other_cofactors
        else:
            raise ValueError("Other Cofactors Component has not been found")

    def set_other_cofactors(self, other_cofactors):
        self._other_cofactors = other_cofactors

    @property
    def has_other_cofactors(self):
        if self._other_cofactors is not None:
            return True
        else:
            return False

    @property
    def get_lipids(self):
        if self._lipids is not None:
            return self._lipids
        else:
            raise ValueError("Lipids Component has not been found")

    def set_lipids(self, lipids):
        self._lipids = lipids

    @property
    def has_lipids(self):
        if self._lipids is not None:
            return True
        else:
            return False

    @property
    def get_nucleics(self):
        if self._nucleics is not None:
            return self._nucleics
        else:
            raise ValueError("Nucleics Component has not been found")

    def set_nucleics(self, nucleics):
        self._nucleics = nucleics

    @property
    def has_nucleids(self):
        if self._nucleics is not None:
            return True
        else:
            return False

    @property
    def get_other_nucleics(self):
        if self._other_nucleics is not None:
            return self._other_nucleics
        else:
            raise ValueError("Other Nucleics Component has not been found")

    def set_other_nucleics(self, other_nucleics):
        self._other_nucleics = other_nucleics

    @property
    def has_other_nucleids(self):
        if self._other_nucleics is not None:
            return True
        else:
            return False

    @property
    def get_components_title(self):
        return self._components_title

    def set_components_title(self, title):
        self._components_title = title

    @property
    def get_components(self):
        return self._components

    def set_component_by_name(self, comp_name, comp):

        if comp_name not in self._component_names:
            raise ValueError("The component name {} is not supported. Allowed: {}".format(comp_name,
                                                                                          self._component_names))
        self._components[comp_name] = comp
