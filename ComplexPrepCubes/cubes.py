from ComplexPrepCubes import utils
from OpenMMCubes import utils as pack_utils
from floe.api import ParallelMixin, parameter
from openeye import oechem
import traceback
from simtk import unit
from simtk.openmm import app
from oeommtools import utils as oeommutils
from oeommtools.packmol import oesolvate
import parmed

from cuberecord import OERecordComputeCube, OEField, OERecord
from cuberecord.ports import RecordInputPort
from datarecord import Types, Meta, ColumnMeta
from cuberecord.constants import DEFAULT_MOL_NAME


class HydrationSetCube(ParallelMixin, OERecordComputeCube):
    title = "Hydration Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem', 'OpenMM', 'PDBFixer']
    description = """
    This cube solvate the molecular system

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the molecular system

    Output:
    -------
    oechem.OEDataRecord - Emits the solvated system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10.0,
        help_text="Padding around protein for solvent box (angstroms)")

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=50.0,
        help_text="Salt concentration (millimolar)")

    ref_structure = parameter.BooleanParameter(
        'ref_structure',
        default=True,
        help_text="If Checked/True the molecule before solvation is attached to the solvated one")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt

            field_system = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                   meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))

            if not record.has_value(field_system):
                self.log.warn("Missing molecule '{}' field".format(field_system.get_name()))
                self.failure.emit(record)
                return

            system = record.get_value(field_system)

            field_system_id = OEField("ID", Types.String)

            if not record.has_value(field_system_id):
                self.log.warn("Missing molecule ID '{}' field".format(field_system_id.get_name()))
                system_id = system.GetTitle()
            else:
                system_id = record.get_value(field_system_id)

            # if opt['ref_structure']:
            #     ref_field_mol = OEField("Reference Molecule", Types.Chem.Mol)
            #     record.set_value(ref_field_mol, system)

            # Solvate the system. Note that the solvated system is translated to the
            # OpenMM cube cell
            sol_system = utils.hydrate(system, opt)
            sol_system.SetTitle(system_id)

            record.set_value(field_system, sol_system)
            record.set_value(field_system_id, system_id)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return


class SolvationSetCube(ParallelMixin, OERecordComputeCube):
    title = "Solvation Cube Packmol"
    version = "0.0.0"
    classification = [["Preparation", "OEChem"]]
    tags = ['OEChem', 'PackMol']
    description = """
    This cube solvate a molecular system

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the molecular system

    Output:
    -------
    oechem.OEDataRecord - Emits the solvated system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    density = parameter.DecimalParameter(
        'density',
        default=1.0,
        help_text="Solution density in g/ml")

    padding_distance = parameter.DecimalParameter(
        'padding_distance',
        default=10.0,
        help_text="The padding distance between the solute and the box edge in A")

    distance_between_atoms = parameter.DecimalParameter(
        'distance_between_atoms',
        default=2.0,
        help_text="The minimum distance between atoms in A")

    solvents = parameter.StringParameter(
        'solvents',
        required=True,
        default='[H]O[H]',
        help_text='Select solvents. The solvents are specified as comma separated smiles strings'
                  'e.g. [H]O[H], C(Cl)(Cl)Cl, CS(=O)C')

    molar_fractions = parameter.StringParameter(
        'molar_fractions',
        default='1.0',
        help_text="Molar fractions of each solvent components. The molar fractions are specified"
                  "as comma separated molar fractions strings e.g. 0.5,0.2,0.3")

    geometry = parameter.StringParameter(
        'geometry',
        default='box',
        choices=['box', 'sphere'],
        help_text="Geometry selection: box or sphere. Sphere cannot be used as periodic system "
                  "along with MD simulation")

    close_solvent = parameter.BooleanParameter(
        'close_solvent',
        default=False,
        help_text="If Checked/True solvent molecules will be placed very close to the solute")

    salt = parameter.StringParameter(
        'salt',
        default='[Na+], [Cl-]',
        help_text='Salt type. The salt is specified as list of smiles strings. '
                  'Each smiles string is the salt component dissociated in the '
                  'solution e.g. Na+, Cl-')

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=0.0,
        help_text="Salt concentration in millimolar")

    neutralize_solute = parameter.BooleanParameter(
        'neutralize_solute',
        default=True,
        help_text='Neutralize the solute by adding Na+ and Cl- counter-ions based on'
                  'the solute formal charge')

    ref_structure = parameter.BooleanParameter(
        'ref_structure',
        default=True,
        help_text="If Checked/True the molecule before solvation is attached to the solvated one")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = self.opt

            field_system = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                   meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))

            if not record.has_value(field_system):
                self.log.warn("Missing molecule '{}' field".format(field_system.get_name()))
                self.failure.emit(record)
                return

            solute = record.get_value(field_system)

            field_system_id = OEField("ID", Types.String)

            if not record.has_value(field_system_id):
                self.log.warn("Missing molecule ID '{}' field".format(field_system_id.get_name()))
                solute_id = solute.GetTitle()
            else:
                solute_id = record.get_value(field_system_id)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(solute) if dp.GetTag() in
                            ["solvents", "molar_fractions", "density"]}
            if new_args:
                for k in new_args:
                    if k == 'molar_fractions':
                        continue
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(solute_id, new_args))
                opt.update(new_args)

            # Solvate the system
            sol_system = oesolvate(solute, **opt)
            self.log.info("Solvated System atom number: {}".format(sol_system.NumAtoms()))
            sol_system.SetTitle(solute_id)

            # if opt['ref_structure']:
            #     ref_field_mol = OEField("Reference Molecule", Types.Chem.Mol)
            #     record.set_value(ref_field_mol, system)

            record.set_value(field_system, sol_system)
            record.set_value(field_system_id, solute_id)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return


class ComplexSetPrepCube(OERecordComputeCube):
    title = "Complex Preparation Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem']
    description = """
        This cube assembles the complex made of the protein and the
        ligands. If a ligand presents multiple conformers, then each conformer
        is bonded to the protein to form the solvated complex. For example if a
        ligand has 3 conformers then 3 complexes are generated.

        Input:
        -------
        oechem.OEDataRecord - Streamed-in of the protein and ligands

        Output:
        -------
        oechem.OEDataRecord - Emits the complex molecules
        """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    remove_explicit_solvent = parameter.BooleanParameter(
        'remove_explicit_solvent',
        default=False,
        description='If True/Checked removes water and ion molecules from the system')

    protein_port = RecordInputPort("protein_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.wait_on('protein_port')
        self.count = 0
        self.check_protein = False

    def process(self, record, port):
        try:
            # protein
            if port == 'protein_port':

                field_mol = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                    meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))

                if not record.has_value(field_mol):
                    self.log.warn("Missing Protein '{}' field".format(field_mol.get_name()))
                    self.failure.emit(record)
                    return

                mol = record.get_value(field_mol)

                field_id = OEField("ID", Types.String)

                if not record.has_value(field_id):
                    self.log.warn("Missing Protein ID '{}' field".format(field_id.get_name()))
                    self.protein_id = 'p' + mol.GetTitle()+'_'+str(self.count)
                else:
                    self.protein_id = record.get_value(field_id)

                # Remove from solution water and ions
                if self.opt['remove_explicit_solvent']:
                    mol = oeommutils.strip_water_ions(mol)

                self.protein = mol
                self.check_protein = True
                return

            # ligands
            if self.check_protein:

                field_mol = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)

                if not record.has_value(field_mol):
                    self.log.warn("Missing Ligand '{}' field".format(field_mol.get_name()))
                    self.failure.emit(record)
                    return

                mol = record.get_value(field_mol)

                field_id = OEField("ID", Types.String)

                if not record.has_value(field_id):
                    self.log.warn("Missing Ligand ID '{}' field".format(field_id.get_name()))
                    lig_id = self.protein_id + '_l' + mol.GetTitle()[0:12] + '_' + str(self.count)
                else:
                    field_id = record.get_value(field_id)
                    lig_id = self.protein_id + '_' + field_id

                num_conf = 0

                for conf in mol.GetConfs():
                    conf_mol = oechem.OEMol(conf)
                    complx = self.protein.CreateCopy()
                    oechem.OEAddMols(complx, conf_mol)

                    # Split the complex in components
                    protein, ligand, water, excipients = oeommutils.split(complx)

                    # If the protein does not contain any atom emit a failure
                    if not protein.NumAtoms():  # Error: protein molecule is empty
                        oechem.OEThrow.Error("The protein molecule does not contains atoms")

                    # If the ligand does not contain any atom emit a failure
                    if not ligand.NumAtoms():  # Error: ligand molecule is empty
                        oechem.OEThrow.Error("The Ligand molecule does not contains atoms")

                    # Check if the ligand is inside the binding site. Cutoff distance 3A
                    if not oeommutils.check_shell(ligand, protein, 3):
                        oechem.OEThrow.Error("The ligand is probably outside the protein binding site")

                    # Removing possible clashes between the ligand and water or excipients
                    if water.NumAtoms():
                        water_del = oeommutils.delete_shell(ligand, water, 1.5, in_out='in')

                    if excipients.NumAtoms():
                        excipient_del = oeommutils.delete_shell(ligand, excipients, 1.5, in_out='in')

                    # Reassemble the complex
                    new_complex = protein.CreateCopy()
                    oechem.OEAddMols(new_complex, ligand)
                    if excipients.NumAtoms():
                        oechem.OEAddMols(new_complex, excipient_del)
                    if water.NumAtoms():
                        oechem.OEAddMols(new_complex, water_del)

                    complex_id_c = lig_id

                    if mol.GetMaxConfIdx() > 1:
                        complex_id_c = complex_id_c + '_c' + str(num_conf)

                    new_complex.SetTitle(complex_id_c)

                    # New Data Record Settings
                    new_complex_field = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                                meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))
                    new_record_complex = OERecord()
                    new_record_complex.set_value(new_complex_field, new_complex)

                    complex_field_id = OEField("ID", Types.String)
                    new_record_complex.set_value(complex_field_id, complex_id_c)

                    num_conf += 1
                    self.success.emit(new_record_complex)
                self.count += 1

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return


class ForceFieldSetCube(ParallelMixin, OERecordComputeCube):
    title = "Force Field Application Cube"
    version = "0.0.0"
    classification = [["Force Field Application", "OEChem"]]
    tags = ['OEChem', 'OEBio', 'OpenMM']
    description = """
    The system is parametrized by using the selected force fields

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the system to parametrize

    Output:
    -------
    oechem.OEDataRecord - Emits the force field parametrized system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_forcefield = parameter.StringParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Force field parameters for protein')

    solvent_forcefield = parameter.StringParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Force field parameters for solvent')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        required=True,
        default='GAFF',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field to parametrize the ligand')

    ligand_res_name = parameter.StringParameter(
        'ligand_res_name',
        required=True,
        default='LIG',
        help_text='Ligand residue name')

    other_forcefield = parameter.StringParameter(
        'other_forcefield',
        required=True,
        default='GAFF',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field used to parametrize other molecules not recognized by the protein force field')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt

            field_system = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)

            if not record.has_value(field_system):
                self.log.warn("Missing molecule '{}' fied".format(field_system.get_name()))
                self.failure.emit(record)
                return

            system = record.get_value(field_system)

            field_system_id = OEField("ID", Types.String)

            if not record.has_value(field_system_id):
                self.log.warn("Missing molecule ID '{}' field".format(field_system_id.get_name()))
                system_id = system.GetTitle()
            else:
                system_id = record.get_value(field_system_id)

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = oeommutils.split(system, ligand_res_name=opt['ligand_res_name'])

            self.log.info("\nComplex name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(system_id,
                                                                                         protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))

            # Unique prefix name used to output parametrization files
            opt['prefix_name'] = system_id

            oe_mol_list = []
            par_mol_list = []

            # Apply FF to the Protein
            if protein.NumAtoms():
                oe_mol_list.append(protein)
                protein_structure = utils.applyffProtein(protein, opt)
                par_mol_list.append(protein_structure)

            # Apply FF to the ligand
            if ligand.NumAtoms():
                oe_mol_list.append(ligand)
                ligand_structure = utils.applyffLigand(ligand, opt)
                par_mol_list.append(ligand_structure)

            # Apply FF to water molecules
            if water.NumAtoms():
                oe_mol_list.append(water)
                water_structure = utils.applyffWater(water, opt)
                par_mol_list.append(water_structure)

            # Apply FF to the excipients
            if excipients.NumAtoms():
                excipient_structure = utils.applyffExcipients(excipients, opt)
                par_mol_list.append(excipient_structure)

                # The excipient order is set equal to the order in related
                # parmed structure to avoid possible atom index mismatching
                excipients = oeommutils.openmmTop_to_oemol(excipient_structure.topology,
                                                           excipient_structure.positions,
                                                           verbose=False)
                oechem.OEPerceiveBondOrders(excipients)
                oe_mol_list.append(excipients)

            # Build the overall Parmed structure
            system_structure = parmed.Structure()

            for struc in par_mol_list:
                system_structure = system_structure + struc

            system_reassembled = oe_mol_list[0].CreateCopy()
            num_atom_system = system_reassembled.NumAtoms()

            for idx in range(1, len(oe_mol_list)):
                oechem.OEAddMols(system_reassembled, oe_mol_list[idx])
                num_atom_system += oe_mol_list[idx].NumAtoms()

            if not num_atom_system == system_structure.topology.getNumAtoms():
                oechem.OEThrow.Error("Parmed and OE topologies mismatch atom number "
                                     "error for system: {}".format(system_id))

            system_reassembled.SetTitle(system_id)

            # Set Parmed structure box_vectors
            is_periodic = True
            try:
                vec_data = pack_utils.PackageOEMol.getData(system_reassembled, tag='box_vectors')
                vec = pack_utils.PackageOEMol.decodePyObj(vec_data)
                system_structure.box_vectors = vec
            except:
                is_periodic = False
                self.log.warn("System {} has been parametrize without periodic box vectors "
                              "for vacuum simulation".format(system_id))

            # Set atom serial numbers, Ligand name and HETATM flag
            for at in system_reassembled.GetAtoms():
                thisRes = oechem.OEAtomGetResidue(at)
                thisRes.SetSerialNumber(at.GetIdx())
                if thisRes.GetName() == 'UNL':
                    # thisRes.SetName("LIG")
                    thisRes.SetHetAtom(True)
                oechem.OEAtomSetResidue(at, thisRes)

            if system_reassembled.GetMaxAtomIdx() != system_structure.topology.getNumAtoms():
                raise ValueError("OEMol system {} and generated Parmed structure "
                                 "mismatch atom numbers".format(system_id))

            # Check if it is possible to create the OpenMM System
            if is_periodic:
                system_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                              nonbondedCutoff=10.0 * unit.angstroms,
                                              constraints=app.HBonds,
                                              removeCMMotion=False)
            else:
                system_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                              constraints=app.HBonds,
                                              removeCMMotion=False)

            record.set_value(field_system, system_reassembled)
            record.set_value(field_system_id, system_id)

            parmed_field = OEField("Parmed", pack_utils.ParmedData)
            record.set_value(parmed_field, system_structure)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return