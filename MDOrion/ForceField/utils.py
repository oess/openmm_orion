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

from openeye import oechem

from MDOrion.ForceField.ff_library import ff_library

from MDOrion.ForceField.ffutils import (ParamMolStructure,
                                        parametrize_component,
                                        clean_tags)

from oeommtools import utils as oeommutils

from oeommtools import data_utils

from simtk.openmm import app

from openeye import oespruce

import parmed

from oemdtoolbox.ForceField import nsr_template_generator

from io import StringIO

from oeommtools.utils import sanitizeOEMolecule

import copy

import numpy as np


class MDComponents:

    def __init__(self, system_representation, components_title="Flask"):

        du = None
        molecules = None

        if type(system_representation) == oechem.OEDesignUnit:
            du = system_representation
        elif type(system_representation) == oechem.OEMol:
            molecules = system_representation
        else:
            raise ValueError("The MDComponent class can be initialized with an OE Design Unit "
                             "or an OE Mol. The object passed is: {}".format(type(system_representation)))

        # MD allowed components
        self._protein = None
        self._ligand = None
        self._other_ligands = None
        self._counter_ions = None
        self._metals = None
        self._excipients = None
        self._solvent = None
        self._water = None
        self._cofactors = None
        self._other_cofactors = None
        self._lipids = None
        self._nucleics = None
        self._other_nucleics = None

        self._component_names = ['protein', 'ligand',
                                 'other_ligands', 'counter_ions',
                                 'metals', 'excipients',
                                 'solvent', 'water', 'cofactors',
                                 'other_cofactors', 'lipids',
                                 'nucleics', 'other_nucleics']

        # Components found in the system_representation
        self._components = dict()

        self._components_title = components_title

        self._box_vectors = None

        if du is not None:
            print("Found DU")
            self._initialize_from_du(du)
        else:
            self._initialize_from_molecules(molecules)

    def __repr__(self):
        ret_str = "\n" + 28 * "-" + "\n"
        ret_str += "{}\n".format(self.get_title)[:28]
        ret_str += "\n{:<20} {:>7}\n".format("Comp_name", "Atoms")
        ret_str += 28 * "-" + "\n"
        for comp_name, comp in self._components.items():
            ret_str += "{:<20} {:>7}\n".format(comp_name, comp.NumAtoms())
        ret_str += 28 * "-" + "\n"
        ret_str += "{:<20} {:>7}\n".format("Total_Atoms", self.num_atoms)

        return ret_str

    def _initialize_from_du(self, du):

        if du.GetTitle() and self._components_title == "Flask":
            self._components_title = du.GetTitle()

        for pair in du.GetTaggedComponents():

            comp_name = pair[0]
            comp_id = du.GetComponentID(comp_name)
            comp = pair[1]

            # Removing Interaction Hint Container and Style from the components
            oechem.OEDeleteInteractionsHintSerializationData(comp)
            oechem.OEDeleteInteractionsHintSerializationIds(comp)
            oechem.OEClearStyle(comp)

            # Clean R-Groups
            for atom in comp.GetAtoms(oechem.OEIsRGroup()):
                nbrs = [nbr for nbr in atom.GetAtoms()]
                comp.DeleteAtom(atom)
                for nbr in nbrs:
                    oechem.OEAssignMDLHydrogens(nbr)
            oechem.OEAddExplicitHydrogens(comp)

            # Reset Comp Order
            oechem.OEPDBOrderAtoms(comp, False)

            if comp_id == oechem.OEDesignUnitComponents_Protein:
                self._protein = comp
                self._components['protein'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Ligand:
                self._ligand = comp
                self._components['ligand'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_OtherLigands:
                self._other_ligands = comp
                self._components['other_ligands'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_CounterIons:
                self._counter_ions = comp
                self._components['counter_ions'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Metals:
                self._metals = comp
                self._components['metals'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Excipients:
                self._excipients = comp
                self._components['excipients'] = comp
            elif comp_id == oechem.OEDesignUnitComponents_Solvent:
                # Separate Water from Spruce Solvent
                pred_water = oechem.OEIsWater(checkHydrogens=True)
                water = oechem.OEMol()
                oechem.OESubsetMol(water, comp, pred_water)
                if water.NumAtoms():
                    self._water = water
                    self._components['water'] = water
                    pred_not_water = oechem.OENotAtom(oechem.OEIsWater(checkHydrogens=True))
                    solvent_not_water = oechem.OEMol()
                    oechem.OESubsetMol(solvent_not_water, comp, pred_not_water)
                    if solvent_not_water.NumAtoms():
                        self._solvent = solvent_not_water
                        self._components['solvent'] = solvent_not_water
                else:
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

        if molecules.GetTitle() and self._components_title == 'Flask':
            self._components_title = molecules.GetTitle()

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
            print("Design Unit Built from Spruce")
        else:
            # Split the complex in components
            protein, ligand, water, excipients = oeommutils.split(molecules,
                                                                  ligand_res_name='LIG')
            if protein.NumAtoms():
                protein = clean_tags(protein)
                oechem.OEPDBOrderAtoms(protein, False)
                self._protein = protein
                self._components['protein'] = protein
            if ligand.NumAtoms():
                ligand = clean_tags(ligand)
                ligand = sanitizeOEMolecule(ligand)
                oechem.OEPDBOrderAtoms(ligand, False)
                self._ligand = ligand
                self._components['ligand'] = ligand
            if water.NumAtoms():
                water = clean_tags(water)
                oechem.OEPDBOrderAtoms(water, False)
                self._water = water
                self._components['water'] = water
            if excipients.NumAtoms():
                excipients = clean_tags(excipients)
                excipients = sanitizeOEMolecule(excipients)
                oechem.OEPDBOrderAtoms(excipients, False)
                self._excipients = excipients
                self._components['excipients'] = excipients

            print("Design Unit Built from Molecule")

        if self.num_atoms != molecules.NumAtoms():
            raise ValueError("Atom number mismatch: {} vs {}".format(self.num_atoms, molecules.NumAstoms()))

    def __getstate__(self):

        def mol_to_bytes(mol):
            return oechem.OEWriteMolToBytes(oechem.OEFormat_OEB, True, mol)

        state = dict(protein=mol_to_bytes(self._protein) if self._protein else None,
                     ligand=mol_to_bytes(self._ligand) if self._ligand else None,
                     other_ligands=mol_to_bytes(self._other_ligands) if self._other_ligands else None,
                     counter_ions=mol_to_bytes(self._counter_ions) if self._counter_ions else None,
                     metals=mol_to_bytes(self._metals) if self._metals else None,
                     excipients=mol_to_bytes(self._excipients) if self._excipients else None,
                     solvent=mol_to_bytes(self._solvent) if self._solvent else None,
                     water=mol_to_bytes(self._water) if self._water else None,
                     cofactors=mol_to_bytes(self._cofactors) if self._cofactors else None,
                     other_cofactors=mol_to_bytes(self._other_cofactors) if self._other_cofactors else None,
                     lipids=mol_to_bytes(self._lipids) if self._lipids else None,
                     nucleics=mol_to_bytes(self._nucleics) if self._nucleics else None,
                     other_nucleics=mol_to_bytes(self._other_nucleics) if self._other_nucleics else None,
                     components_title=self._components_title,
                     box_vectors=data_utils.encodePyObj(self._box_vectors) if self._box_vectors else None
                     )

        return state

    def __setstate__(self, state):

        # Supported component Names
        self._component_names = ['protein', 'ligand',
                                 'other_ligands', 'counter_ions',
                                 'metals', 'excipients',
                                 'solvent', 'water', 'cofactors',
                                 'other_cofactors', 'lipids',
                                 'nucleics', 'other_nucleics']

        def mol_from_bytes(mol_bytes):
            mol = oechem.OEMol()
            oechem.OEReadMolFromBytes(mol, oechem.OEFormat_OEB, True, mol_bytes)
            return mol

        self._components = dict()
        # self._total_atoms = 0

        for comp_name, comp in state.items():

            if comp_name == 'components_title':
                self._components_title = comp
                continue

            if comp_name == 'box_vectors':
                if comp is not None:
                    box_vec = data_utils.decodePyObj(comp)
                    self._box_vectors = box_vec
                else:
                    self._box_vectors = None
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
            elif comp_name == 'water':
                self._water = mol
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
    def copy(self):
        data = copy.deepcopy(self)
        return data

    @property
    def num_atoms(self):

        tot_atoms = 0
        for comp_name, comp in self._components.items():
            tot_atoms += comp.NumAtoms()

        return tot_atoms

    @property
    def create_flask(self):
        flask = oechem.OEMol()
        map_dic = dict()
        for comp_name, comp in self._components.items():

            # amap is a list. The index i in the list is the atom Idx in the
            # source molecule and the corresponding list element at position i
            # is the OE Atom in the destination molecule
            amap, bmap = oechem.OEAddMols(flask, comp)

            if not amap:
                raise ValueError("The flask cannot be created. Problems with the component: {}".format(comp_name))

            map_dic[comp_name] = [at.GetIdx() for at in amap]

        flask.SetTitle(self._components_title)
        return flask, map_dic

    def update_components_coords(self, flask_new_coords):

        flask, map_dic = self.create_flask

        for at_flask, at_flask_nc in zip(flask.GetAtoms(), flask_new_coords.GetAtoms()):
            if at_flask.GetName() != at_flask_nc.GetName():
                raise ValueError("The provided topology in not in sync with the MD Components topology")

        new_coords = flask_new_coords.GetCoords()

        for comp_name, oe_comp in self.get_components.items():

            with oechem.oemolostream(comp_name+'.oeb') as ofs:
                oechem.OEWriteConstMolecule(ofs, oe_comp)

            idx_comp = map_dic[comp_name]
            coords = []
            for at in oe_comp.GetAtoms():
                coords.append(new_coords[idx_comp[at.GetIdx()]])

            oe_comp.SetCoords(oechem.OEFloatArray(np.array(coords).ravel()))

            self.set_component_by_name(comp_name, oe_comp)

            # if comp_name == 'protein':
            #     self.set_protein(oe_comp)
            # elif comp_name == 'ligand':
            #     self.set_ligand(oe_comp)
            # elif comp_name == 'other_ligands':
            #     self.set_other_ligands(oe_comp)
            # elif comp_name == 'counter_ions':
            #     self.set_counter_ions(oe_comp)
            # elif comp_name == 'metals':
            #     self.set_metals(oe_comp)
            # elif comp_name == 'excipients':
            #     self.set_excipients(oe_comp)
            # elif comp_name == 'solvent':
            #     self.set_solvent(oe_comp)
            # elif comp_name == 'water':
            #     self.set_water(oe_comp)
            # elif comp_name == 'cofactors':
            #     self.set_cofactors(oe_comp)
            # elif comp_name == 'other_cofactors':
            #     self.set_other_cofactors(oe_comp)
            # elif comp_name == 'lipids':
            #     self.set_lipids(oe_comp)
            # elif comp_name == 'nucleics':
            #     self.set_nucleics(oe_comp)
            # elif comp_name == 'other_nucleics':
            #     self.set_other_nucleics(oe_comp)

    @property
    def get_protein(self):
        if self._protein is not None:
            return self._protein
        else:
            raise ValueError("Protein Component has not been found")

    def set_protein(self, protein):
        self._protein = protein
        self._components['protein'] = self._protein

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
        self._components['ligand'] = self._ligand

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
        self._components['other_ligands'] = self._other_ligands

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
        self._components['counter_ions'] = self._counter_ions

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
        self._components['metals'] = self._metals

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
        self._components['excipients'] = self._excipients

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
        self._components['solvent'] = self._solvent

    @property
    def has_solvent(self):
        if self._solvent is not None:
            return True
        else:
            return False

    @property
    def get_water(self):
        if self._water is not None:
            return self._water
        else:
            raise ValueError("Water Component has not been found")

    def set_water(self, water):
        self._water = water
        self._components['water'] = self._water

    @property
    def has_water(self):
        if self._water is not None:
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
        self._components['cofactors'] = self._cofactors

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
        self._components['other_cofactors'] = self._other_cofactors

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
        self._components['lipids'] = self._lipids

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
        self._components['nucleics'] = self._nucleics

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
        self._components['other_nucleics'] = self._other_nucleics

    @property
    def has_other_nucleids(self):
        if self._other_nucleics is not None:
            return True
        else:
            return False

    @property
    def has_box_vectors(self):
        if self._box_vectors is not None:
            return True
        else:
            return False

    @property
    def get_box_vectors(self):
        if self._box_vectors is not None:
            return self._box_vectors
        else:
            raise ValueError("Box Vectors has not been found")

    def set_box_vectors(self, box_vectors):
        self._box_vectors = box_vectors

    @property
    def get_title(self):
        return self._components_title

    def set_title(self, title):
        self._components_title = title

    @property
    def get_info(self):
        return self.__repr__()

    @property
    def get_components(self):
        return self._components

    def set_component_by_name(self, comp_name, comp):

        if comp_name not in self._component_names:
            raise ValueError("The component name {} is not supported. Allowed: {}".format(comp_name,
                                                                                          self._component_names))
        self._components[comp_name] = comp

        if comp_name == 'protein':
            self._protein = comp
        elif comp_name == 'ligand':
            self._ligand = comp
        elif comp_name == 'other_ligands':
            self._other_ligands = comp
        elif comp_name == 'counter_ions':
            self._counter_ions = comp
        elif comp_name == 'metals':
            self._metals = comp
        elif comp_name == 'excipients':
            self._excipients = comp
        elif comp_name == 'solvent':
            self._solvent = comp
        elif comp_name == 'water':
            self._water = comp
        elif comp_name == 'cofactors':
            self._cofactors = comp
        elif comp_name == 'other_cofactors':
            self._other_cofactors = comp
        elif comp_name == 'lipids':
            self._lipids = comp
        elif comp_name == 'nucleics':
            self._nucleics = comp
        elif comp_name == 'other_nucleics':
            self._other_nucleics = comp

    def parametrize_components(self, protein_ff='Amber14SB',
                               ligand_ff='OpenFF_1.0.0',
                               other_ff='OpenFF_1.0.0'):

        ParamMDComp = ParametrizeMDComponents(self, protein_ff, ligand_ff, other_ff)
        pmd = ParamMDComp.parametrize_components

        if self.has_box_vectors:
            pmd.box_vectors = self.get_box_vectors

        self.ParamMDComp = ParamMDComp

        return pmd


class ParametrizeMDComponents:
    def __init__(self, md_component,
                 protein_ff='Amber14SB',
                 ligand_ff='OpenFF_1.2.0',
                 other_ff='OpenFF_1.2.0'):

        # MD Components
        self.md_components = md_component

        # Force Field to use: User Defined
        self.protein_ff = ff_library.proteinff[protein_ff]
        self.ligand_ff = ff_library.ligandff[ligand_ff]
        self.other_ff = ff_library.otherff[other_ff]

        # Extended Force Field
        # self.protein_extended_ff = ff_library.protein_extended_ff[protein_ff]

        # Force Field to use: Default
        self.counter_ions_ff = ff_library.counter_ionsff['Counter_ions']
        self.metals_ff = ff_library.metals_ff['Metals']
        self.excipients_ff = ff_library.excipients_ff['Excipients']
        self.solvent_ff = ff_library.solventff['Tip3p']
        self.cofactors_ff = ff_library.cofactors_ff['Cofactors']
        self.lipids_ff = ff_library.lipids_ff['Lipids']
        self.nucleics_ff = ff_library.nucleics_ff['Nucleics']

    @staticmethod
    def _check_formal_vs_partial_charge(comp_name, component, pmd_component):
        formal_charge = 0
        for at in component.GetAtoms():
            formal_charge += at.GetFormalCharge()

        partial_charge = 0.0
        for at in pmd_component.atoms:
            partial_charge += at.charge

        if abs(formal_charge - partial_charge) > 0.01:
            raise ValueError("Component: {} - Formal charge and Parmed Partial charge mismatch: {} vs {}".
                             format(comp_name, formal_charge, partial_charge))

    @property
    def parametrize_protein(self):
        if self.md_components.has_protein:
            print("Protein Parametrized by using the ff: {}".format(self.protein_ff))
            protein = self.md_components.get_protein

            # OpenMM topology and positions from OEMol
            topology, positions = oeommutils.oemol_to_openmmTop(protein)

            # Try to apply the selected FF on the Protein
            forcefield = app.ForceField(self.protein_ff)

            unmatched_res_list = forcefield.getUnmatchedResidues(topology)

            # If there are force field unmatched residues
            if unmatched_res_list:

                standard_resides_fail = []

                for res in unmatched_res_list:
                    if res.name in ff_library.protein_standard_residue_names:
                        standard_resides_fail.append(res)

                if standard_resides_fail:

                    omm_residues = []
                    for res in topology.residues():
                        omm_residues.append(res)

                    oe_residues = []
                    for oe_res in oechem.OEGetResidues(protein):
                        oe_residues.append(oe_res)

                    map_omm_oe = dict()
                    for omm_res, oe_res in zip(omm_residues, oe_residues):
                        map_omm_oe[omm_res] = oe_res

                    info_fail = [(map_omm_oe[res].GetName(),
                                  map_omm_oe[res].GetResidueNumber(),
                                  map_omm_oe[res].GetChainID())
                                 for res in standard_resides_fail]

                    raise ValueError(
                        "The following protein residues cannot be parametrize: {}".format(info_fail))

                # Try to parametrize Non Standard Residues
                ffxml_nsr_template_list = nsr_template_generator(protein,
                                                                 topology,
                                                                 forcefield,
                                                                 openff=ff_library.ligandff['OpenFF_1.0.0'])

                for ffxml_template in ffxml_nsr_template_list:
                    forcefield.loadFile(StringIO(ffxml_template))

            omm_protein = forcefield.createSystem(topology, rigidWater=False, constraints=None)
            protein_pmd = parmed.openmm.load_topology(topology, omm_protein, xyz=positions)

            self._check_formal_vs_partial_charge("protein", self.md_components.get_protein, protein_pmd)

            return protein_pmd
        else:
            raise ValueError("Protein is not present in the MDComponents")

    @property
    def parametrize_ligand(self):

        if self.md_components.has_ligand:
            print("Ligand Parametrized by using the ff: {}".format(self.ligand_ff))
            prefix_name = 'LIG'

            pmd = ParamMolStructure(self.md_components.get_ligand,
                                    self.ligand_ff,
                                    prefix_name=prefix_name,
                                    recharge=False)

            ligand_pmd = pmd.parameterize()
            ligand_pmd.residues[0].name = prefix_name

            self._check_formal_vs_partial_charge("ligand", self.md_components.get_ligand, ligand_pmd)

            return ligand_pmd
        else:
            raise ValueError("Ligand is not present in the MDComponents")

    @property
    def parametrize_other_ligands(self):
        print("Other Ligands Parametrized by using the ff: {}".format(self.ligand_ff))
        if self.md_components.has_other_ligands:

            other_ligand_pmd = parametrize_component(self.md_components.get_other_ligands,
                                                     self.protein_ff,
                                                     self.ligand_ff)

            self._check_formal_vs_partial_charge("other ligands", self.md_components.get_other_ligands, other_ligand_pmd)

            return other_ligand_pmd
        else:
            raise ValueError("Other Ligands are not present in the MDComponents")

    @property
    def parametrize_counter_ions(self):

        if self.md_components.has_counter_ions:
            counter_ions_pmd = parametrize_component(self.md_components.get_counter_ions,
                                                     self.counter_ions_ff,
                                                     self.other_ff)
            return counter_ions_pmd
        else:
            raise ValueError("Counter Ions are not present in the MDComponents")

    @property
    def parametrize_metals(self):

        if self.md_components.has_metals:
            metals_pmd = parametrize_component(self.md_components.get_metals,
                                               self.metals_ff,
                                               self.other_ff)

            self._check_formal_vs_partial_charge("metals", self.md_components.get_metals, metals_pmd)

            return metals_pmd
        else:
            raise ValueError("Metals are not present in the MDComponents")

    @property
    def parametrize_excipients(self):

        if self.md_components.has_excipients:
            excipients_pmd = parametrize_component(self.md_components.get_excipients,
                                                   self.excipients_ff,
                                                   self.other_ff)

            self._check_formal_vs_partial_charge("excipients", self.md_components.get_excipients, excipients_pmd)

            return excipients_pmd

        else:
            raise ValueError("Excipients are not present in the MDComponents")

    @property
    def parametrize_solvent(self):

        if self.md_components.has_solvent:

            solvent_pmd = parametrize_component(self.md_components.get_solvent,
                                                self.solvent_ff,
                                                self.other_ff)

            self._check_formal_vs_partial_charge("solvent", self.md_components.get_solvent, solvent_pmd)

            return solvent_pmd

        else:
            raise ValueError("Solvent is not present in the MDComponents")

    @property
    def parametrize_water(self):
        if self.md_components.has_water:

            # OpenMM topology and positions from OEMol
            topology, positions = oeommutils.oemol_to_openmmTop(self.md_components.get_water)

            # Try to apply the selected FF on the component
            forcefield = app.ForceField(self.solvent_ff)

            # List of the unrecognized component
            unmatched_res_list = forcefield.getUnmatchedResidues(topology)

            if not unmatched_res_list:
                water_omm_system = forcefield.createSystem(topology, rigidWater=False, constraints=None)
                water_pmd = parmed.openmm.load_topology(topology, water_omm_system, xyz=positions)

                self._check_formal_vs_partial_charge("water", self.md_components.get_water, water_pmd)

                return water_pmd

            else:
                raise ValueError("Water cannot be parametrized by using the FF: {}\n Problematic Residues are: {}".
                                 format(self.solvent_ff, unmatched_res_list))
        else:
            raise ValueError("Water is not present in the MDComponents")

    @property
    def parametrize_cofactors(self):

        if self.md_components.has_cofactors:

            ff = [self.cofactors_ff, self.metals_ff]

            cofactors_pmd = parametrize_component(self.md_components.get_cofactors,
                                                  ff,
                                                  self.other_ff)

            self._check_formal_vs_partial_charge("cofactors", self.md_components.get_cofactors, cofactors_pmd)

            return cofactors_pmd
        else:
            raise ValueError("Cofactors are not present in the MDComponents")

    @property
    def parametrize_other_cofactors(self):

        if self.md_components.has_other_cofactors:

            ff = [self.cofactors_ff, self.metals_ff]

            other_cofactors_pmd = parametrize_component(self.md_components.get_other_cofactors,
                                                        ff,
                                                        self.other_ff)

            self._check_formal_vs_partial_charge("other cofactors", self.md_components.get_other_cofactors,
                                                 other_cofactors_pmd)

            return other_cofactors_pmd
        else:
            raise ValueError("Other Cofactors not present in the MDComponents")

    @property
    def parametrize_lipids(self):

        if self.md_components.has_lipinds:
            lipids_pmd = parametrize_component(self.md_components.get_lipids,
                                               self.lipids_ff,
                                               self.other_ff)

            self._check_formal_vs_partial_charge("lipids", self.md_components.get_lipids, lipids_pmd)

            return lipids_pmd
        else:
            raise ValueError("Lipids are not present in the MDComponents")

    @property
    def parametrize_nucleics(self):

        if self.md_components.has_nucleics:
            nucleics_pmd = parametrize_component(self.md_components.get_nucleics,
                                                 self.nucleics_ff,
                                                 self.other_ff)

            self._check_formal_vs_partial_charge("lipids", self.md_components.get_nucleics, nucleics_pmd)

            return nucleics_pmd
        else:
            raise ValueError("Nucleics molecule not present in the MDComponents")

    @property
    def parametrize_other_nucleics(self):

        if self.md_components.has_other_nucleics:
            other_nucleics_pmd = parametrize_component(self.md_components.get_other_nucleics,
                                                       self.nucleics_ff,
                                                       self.other_ff)
            return other_nucleics_pmd
        else:
            raise ValueError("Other Nucleics molecules are not present in the MDComponents")

    @property
    def parametrize_components(self):

        flask_pmd = parmed.Structure()

        for comp_name, comp in self.md_components.get_components.items():

            if comp_name == 'protein':
                flask_pmd += self.parametrize_protein
            elif comp_name == 'ligand':
                flask_pmd += self.parametrize_ligand
            elif comp_name == 'other_ligands':
                flask_pmd += self.parametrize_other_ligands
            elif comp_name == 'counter_ions':
                flask_pmd += self.parametrize_counter_ions
            elif comp_name == 'metals':
                flask_pmd += self.parametrize_counter_ions
            elif comp_name == 'excipients':
                flask_pmd += self.parametrize_excipients
            elif comp_name == 'solvent':
                flask_pmd += self.parametrize_solvent
            elif comp_name == 'water':
                flask_pmd += self.parametrize_water
            elif comp_name == 'cofactors':
                flask_pmd += self.parametrize_cofactors
            elif comp_name == 'other_cofactors':
                flask_pmd += self.parametrize_other_cofactors
            elif comp_name == 'lipids':
                flask_pmd += self.parametrize_lipids
            elif comp_name == 'nucleics':
                flask_pmd += self.parametrize_nucleics
            elif comp_name == 'other_nucleics':
                flask_pmd += self.parametrize_other_nucleics
            else:
                raise ValueError("The parametrization of the component {} is not supported".format(comp_name))

        return flask_pmd
