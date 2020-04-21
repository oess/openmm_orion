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
                                        parametrize_unknown_component)
import parmed


class ParametrizeDU:
    def __init__(self, design_unit, opt):

        # Design Unit
        self.du = design_unit
        self.opt = opt

        # Force Field to use: User Defined
        self.other_ff = ff_library.otherff[opt['other_forcefield']]
        self.ligand_ff = ff_library.ligandff[opt['ligand_forcefield']]
        self.protein_ff = ff_library.proteinff[opt['protein_forcefield']]

        # Force Field to use: Default
        self.counter_ions_ff = ff_library.counter_ionsff['Counter_ions']
        self.metals_ff = ff_library.metals_ff['Metals']
        self.excipients_ff = ff_library.excipients_ff['Excipients']
        self.solvent_ff = ff_library.solventff['Tip3p']
        self.cofactors_ff = ff_library.cofactors_ff['Cofactors']
        self.lipids_ff = ff_library.lipids_ff['Lipids']
        self.nucleics_ff = ff_library['Nucleics']

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
                opt['Logger'].warn("The component {} parametrization is not currently supported".
                                   format(comp_name))
        if not self.components:
            raise ValueError("None of the DU components cannot be parametrized")

    def parametrize_protein(self):
        # if self.protein:
        #     pmd_protein, unrec_prot = parametrize_component(self.protein, self.protein_ff)
        # else:
        #     raise ValueError("Protein is not present in the DU")
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

                    # self.opt['Logger'].warn("Some Counter Ions have  been parametrized by using the ff: {}".
                    #                         format(self.other_ff))

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

                    # self.opt['Logger'].warn("Some Metal Ions have been parametrized by using the ff: {}".
                    #                         format(self.other_ff))
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

                    # self.opt['Logger'].warn("Some Excipients have been parametrized by using the ff: {}".
                    #                         format(self.other_ff))
                return excipients_pmd
            else:
                raise ValueError("Excipients cannot be parametrized by using the FF: {}".format(self.excipients_ff))
        else:
            raise ValueError("Excipients are not present in the DU")

    @property
    def parametrize_solvent(self):

        if self.solvent:

            # PDB_Order Error. Call Perceive Residue for the wrong PDB Atom Order
            # for oe_res in oechem.OEGetResidues(component_copy):
            #     print(oe_res.GetName(), oe_res.GetResidueNumber())
            #     for oe_at in oechem.OEGetResidueAtoms(component, oe_res, oechem.OEAssumption_BondedResidue +
            #                         oechem.OEAssumption_ResPerceived +
            #                        oechem.OEAssumption_PDBOrder):
            #         print("\t {}".format(oe_at))

            oechem.OEPerceiveResidues(self.solvent, oechem.OEPreserveResInfo_All)

            solvent_pmd, unrec_solvent = parametrize_component(self.solvent, self.solvent_ff)

            if solvent_pmd is not None:

                if unrec_solvent is not None:

                    unk_solvent_pmd = parametrize_unknown_component(unrec_solvent, self.other_ff)
                    solvent_pmd += unk_solvent_pmd

                # self.opt['Logger'].warn("Some Solvent molecules have been parametrized by using the ff: {}".
                #                         format(self.other_ff))
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
                continue
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

