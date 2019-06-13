# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
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

import os
import tempfile

from floe.api import (SourceCube,
                      ParallelMixin, parameter)
from floe.api.orion import in_orion

from datarecord import (OERecord,
                        OEWriteRecord,
                        OEPrimaryMolField)
from datarecord.utils import TemporaryPath

from floe.api.cubes import (ComputeCube,
                            SinkCube)

from orionplatform.mixins import RecordPortsMixin
from orionplatform.ports import RecordOutputPort

from orionclient.session import APISession
from orionclient.types import (Dataset,
                               File)

from openeye import oechem, oedocking

from spruce.prep.options import (OESpruceSplitOptions,
                                 OESprucePrepOptions)
from spruce.prep.iridium_score import (OECalculateDPI,
                                       OECalcIridiumData)
from spruce.pdb.metadata import OECreateStructureMetadata
from spruce.utils.file_utils import sanitize_filename

from .fields import SpruceFields
from .params import SpruceParameters

from .prep_utils import (create_du_ird,
                         DesignUnitToReceptor,
                         download_pdb,
                         download_mtz,
                         get_bu_ref_mol,
                         get_du_ref_file,
                         get_orion_job_tag,
                         get_mtz_url,
                         get_pdb_url,
                         ird_to_str,
                         ReadProteinFromPDB,
                         write_mol_to_bytes,
                         write_mtz)


class PDBCodeToUrl(SourceCube):
    title = "PDB code to URL"
    classification = [["Spruce", "PDB code to URL Translation"]]
    tags = ["OpenEye", "snowball", "Spruce"]
    description = "PDB/MTZ urls are constructed from an input PDB code, or a text file containing multiple codes."  # noqa

    pdb_code = SpruceParameters().pdb_code
    pdb_codes = SpruceParameters().pdb_codes
    design_ref = SpruceParameters().design_ref
    biounit_ref = SpruceParameters().biounit_ref

    success = RecordOutputPort('success')
    failure = RecordOutputPort('failure')

    def _set_record(self, pdb_code, pdb_url, mtz_url):
        record = OERecord()
        record.set_value(SpruceFields().pdb_code, pdb_code)
        record.set_value(SpruceFields().pdb_url, pdb_url)
        record.add_field(SpruceFields().mtz_url)
        if mtz_url:
            record.set_value(SpruceFields().mtz_url, mtz_url)

        record.add_field(SpruceFields().du_ref_file)
        if self.args.design_ref:
            design_ref = self.args.design_ref
            du_ref = get_du_ref_file(design_ref)
            if du_ref:
                du_ref_bytes = oechem.OEWriteDesignUnitToBytes(du_ref)
                record.set_value(SpruceFields().du_ref_file, du_ref_bytes)
            else:
                msg = "Could not read the DU ref file: {}. Aborting.".format(self.args.design_ref)
                return record, False, msg

        record.add_field(SpruceFields().bu_ref_mol)
        if self.args.biounit_ref:
            bu_ref = get_bu_ref_mol(self.args.biounit_ref)
            if bu_ref:
                record.set_value(SpruceFields().bu_ref_mol, bu_ref)
            else:
                msg = " Could not download the BU ref: {}".format(self.args.biounit_ref)
                return record, False, msg

        return record, True, ""

    def __iter__(self):
        if not self.args.pdb_codes and not self.args.pdb_code:
            self.log.info(" A pdb code or a file of PDB codes is required.")  # noqa
            self.failure.emit(OERecord())
            return

        if self.args.pdb_codes and self.args.pdb_code:
            self.log.info(" Only a single PDB code or a file of PDB codes is required.")  # noqa
            self.failure.emit(OERecord())
            return

        if self.args.pdb_codes:
            if not os.path.isfile(str(self.args.pdb_codes)):
                if not isinstance(self.args.pdb_codes, dict):
                    file = self.args.pdb_codes
                else:
                    file = self.args.pdb_codes["file"]
                resource = APISession.get_resource(File, file)
                self.temp = TemporaryPath(suffix=resource.name)
                resource.download_to_file(self.temp.path)
                self.args.pdb_codes = self.temp.path
            try:
                with open(self.args.pdb_codes) as pdb_codes:
                    for pdb_code_raw in pdb_codes.readlines():
                        pdb_code = pdb_code_raw.strip()
                        pdb_code = pdb_code.replace(".pdb", "")
                        if len(pdb_code) == 4:
                            pdb_url = get_pdb_url(pdb_code)
                            mtz_url = get_mtz_url(pdb_code, self.log.info)
                            record, good_record, msg = self._set_record(pdb_code, pdb_url, mtz_url)
                            if not good_record:
                                self.log.info(msg)
                                self.failure.emit(record)
                                return
                            else:
                                yield record
                        else:
                            self.log.info("Skipping this invalid code in the pdb codes file: {}".format(pdb_code))
                            self.failure.emit(OERecord())
            except FileNotFoundError:
                self.log.info(" The pdb code file {} could not be read.".format(self.args.pdb_codes))  # noqa
                self.failure.emit(OERecord())
                return

        else:
            pdb_url = get_pdb_url(self.args.pdb_code)
            mtz_url = get_mtz_url(self.args.pdb_code, self.log.info)

            record, good_record, msg = self._set_record(self.args.pdb_code, pdb_url, mtz_url)
            if not good_record:
                self.log.info(msg)
                self.failure.emit(record)
                return
            else:
                yield record


class UrlToFile(RecordPortsMixin, ComputeCube):
    title = "PDB/MTZ Url to file(s) on an OERecord"
    classification = [["Spruce", "Utilities", "URL to File Translation"]]
    tags = ["OpenEye", "snowball", "Spruce"]
    description = "PDB/MTZ files are downloaded from urls coming from input datarecords."  # noqa

    def process(self, record, port):
        pdb_code = record.get_value(SpruceFields().pdb_code)
        pdb_url = record.get_value(SpruceFields().pdb_url)
        try:
            pdb = download_pdb(pdb_url)
        except IOError:
            self.log.info("PDB ({}) Could not download the PDB. Aborting.".format(pdb_code))  # noqa
            self.failure.emit(record)
            return

        out_record = OERecord()
        meta_data = str(OECreateStructureMetadata(pdb, from_stream=True))
        out_record.add_field(SpruceFields().meta_data)
        if meta_data:
            out_record.set_value(SpruceFields().meta_data, meta_data)

        mtz_url = None
        if record.has_value(SpruceFields().mtz_url):
            mtz_url = record.get_value(SpruceFields().mtz_url)
        mtz = download_mtz(mtz_url, pdb_code, self.log.info)
        out_record.add_field(SpruceFields().mtz_file)
        if mtz:
            out_record.set_value(SpruceFields().mtz_file, mtz)

        out_record.set_value(SpruceFields().pdb_code, pdb_code)
        out_record.set_value(SpruceFields().pdb_file, pdb)

        out_record.add_field(SpruceFields().du_ref_file)
        if record.has_value(SpruceFields().du_ref_file):
            du_ref_bytes = record.get_value(SpruceFields().du_ref_file)
            out_record.set_value(SpruceFields().du_ref_file, du_ref_bytes)

        out_record.add_field(SpruceFields().bu_ref_mol)
        if record.has_value(SpruceFields().bu_ref_mol):
            bu_ref = record.get_value(SpruceFields().bu_ref_mol)
            out_record.set_value(SpruceFields().bu_ref_mol, bu_ref)

        self.success.emit(out_record)


class DUtoReceptorDataset(RecordPortsMixin, SinkCube):
    title = "DU to Receptor Dataset"
    classification = [["Spruce", "DU to Receptor Dataset"]]
    tags = ["OpenEye", "snowball", "Spruce", "I/O"]
    categories = [["I/O", "Writers"]]
    description = "Datasets are created from an input vector of OEDesignUnits on a datarecord."  # noqa

    cache_grids = SpruceParameters().cache_grids
    keep_du = SpruceParameters().keep_du
    dataset_prefix = SpruceParameters().dataset_prefix
    dataset_name = SpruceParameters().dataset_name
    primary_mol_field_name = SpruceParameters().primary_mol_field_name

    def begin(self):
        if in_orion():
            self.in_orion = True
            self.job_tag = get_orion_job_tag()
        else:
            self.in_orion = False
            self.db_file = self.args.dataset_prefix + "{}.oedb"

    def write(self, record, port):
        du_byte_list = record.get_value(SpruceFields().du_vec)

        self.pdb_code = record.get_value(SpruceFields().pdb_code)

        for du_bytes in du_byte_list:
            du = oechem.OEDesignUnit()
            oechem.OEReadDesignUnitFromBytes(du, du_bytes)

            receptor = DesignUnitToReceptor(du, cache_grids=self.args.cache_grids)  # noqa
            if not receptor:
                self.log.warn("Could not make a receptor for DU {}".format(du.GetTitle()))  # noqa
                continue

            out_record = OERecord()
            primary_mol_field = OEPrimaryMolField(default_name=self.args.primary_mol_field_name)
            if self.args.primary_mol_field_name == SpruceFields().ligand.get_name():
                out_record.set_value(primary_mol_field, du.GetImpl().GetLigand())
                out_record.set_value(SpruceFields().receptors, receptor)
            else:
                out_record.set_value(primary_mol_field, receptor)
                out_record.set_value(SpruceFields().ligand, du.GetImpl().GetLigand())

            out_record.set_value(SpruceFields().du_single, du_bytes)
            du_title = du.GetTitle()
            pdb_code = du_title[:du_title.find("(")] if du_title.find("(") != -1 else "Unknown"  # noqa
            out_record.set_value(SpruceFields().pdb_code, pdb_code)
            sq = du.GetStructureQuality()
            iridium_score = oechem.OEGetIridiumCategoryName(sq.GetIridiumData().GetCategory()) if sq.HasIridiumData() else "N/A"  # noqa
            out_record.set_value(SpruceFields().ird_single, iridium_score)
            out_record.set_value(SpruceFields().du_title, du_title)

            fname = 'Receptor_' + self.pdb_code

            if self.in_orion:
                dataset = Dataset.create(APISession, fname)  # noqa
                dataset.write(out_record)
                dataset.finalize()
                APISession.tag_resource(dataset, "Spruce")
                APISession.tag_resource(dataset, "OEDesignUnit")
                APISession.tag_resource(dataset, "Receptor")
                APISession.tag_resource(dataset, du_title)
                APISession.tag_resource(dataset, "PDB: {}".format(pdb_code))
                APISession.tag_resource(dataset, "Iridium: {}".format(iridium_score))
                APISession.tag_resource(dataset, self.job_tag)
            else:

                fname += '.oedb'

                record_stream = oechem.oeofstream(fname)
                OEWriteRecord(record_stream, out_record, fmt='binary')
                record_stream.close()


class DUtoMDDataset(RecordPortsMixin, SinkCube):
    title = "DU to MD Dataset"
    classification = [["Spruce", "DU to MD Dataset"]]
    tags = ["OpenEye", "snowball", "Spruce", "I/O"]
    categories = [["I/O", "Writers", "MD"]]
    description = "MD Datasets are created from an input vector of OEDesignUnits on a datarecord."  # noqa

    success = RecordOutputPort('success')
    failure = RecordOutputPort('failure')

    keep_du = SpruceParameters().keep_du  # TODO: do something with this param

    def begin(self):
        if in_orion():
            self.in_orion = True
            self.job_tag = get_orion_job_tag()
        else:
            self.in_orion = False

        self.component_mask = oechem.OEDesignUnitComponents_Default ^ oechem.OEDesignUnitComponents_Metals ^ oechem.OEDesignUnitComponents_Ligand

    def write(self, record, port):
        du_byte_list = record.get_value(SpruceFields().du_vec)

        self.pdb_code = record.get_value(SpruceFields().pdb_code)

        for du_bytes in du_byte_list:
            du = oechem.OEDesignUnit()
            oechem.OEReadDesignUnitFromBytes(du, du_bytes)
            md_complex = oechem.OEGraphMol()
            if not du.GetComponents(md_complex, self.component_mask):
                self.log.warn("Could not make an MD ready complex from DU named {}".format(du.GetTitle()))  # noqa
                continue

            md_complex.SetTitle(self.pdb_code)

            out_record = OERecord()
            out_record.set_value(SpruceFields().md_complex, md_complex)

            # # TODO DEBUG
            # with oechem.oemolostream("test.oeb") as ofs:
            #     oechem.OEWriteConstMolecule(ofs, md_complex)

            du_title = du.GetTitle()

            pdb_code = du_title[:du_title.find("(")] if du_title.find("(") != -1 else "Unknown"  # noqa
            sq = du.GetStructureQuality()
            iridium_score = oechem.OEGetIridiumCategoryName(sq.GetIridiumData().GetCategory()) if sq.HasIridiumData() else "N/A"  # noqa

            if "biounit" in du_title:
                fname = "BU_" + self.pdb_code

            else:
                fname = "MDReady_" + self.pdb_code

            if self.in_orion:

                dataset = Dataset.create(APISession, fname)  # noqa
                dataset.write(out_record)
                dataset.finalize()
                APISession.tag_resource(dataset, "Spruce")
                APISession.tag_resource(dataset, "OEDesignUnit")
                APISession.tag_resource(dataset, "MD Ready")
                APISession.tag_resource(dataset, du_title)
                APISession.tag_resource(dataset, "PDB: {}".format(pdb_code))
                APISession.tag_resource(dataset, "Iridium: {}".format(iridium_score))
                APISession.tag_resource(dataset, self.job_tag)
            else:
               
                fname += '.oedb'

                record_stream = oechem.oeofstream(fname)
                OEWriteRecord(record_stream, out_record, fmt='binary')
                record_stream.close()


class DUVectorToSingleDU(RecordPortsMixin, ComputeCube):
    title = "A DU vector filters down to a single DU"
    classification = [["Spruce", "DU vector to single DU"]]
    tags = ["OpenEye", "snowball", "Spruce"]
    description = "A single OEDesignUnit is chosen from an input vector of OEDesignUnits."  # noqa

    wanted_du = SpruceParameters().wanted_du

    def process(self, record, port):
        du_byte_list = record.get_value(SpruceFields().du_vec)
        chosen_du = None
        for du_bytes in du_byte_list:
            du = oechem.OEDesignUnit()
            oechem.OEReadDesignUnitFromBytes(du, du_bytes)
            if self.args.wanted_du and self.args.wanted_du in du.GetTitle():
                chosen_du = du_bytes
                break
            else:
                sq = du.GetStructureQuality()
                ird = sq.GetIridiumData()
                if sq.HasIridiumData() and ird.GetCategory() == oechem.OEIridiumCategory_HT:
                    chosen_du = du_bytes
                    break
                elif chosen_du is None and sq.HasIridiumData() and ird.GetCategory() == oechem.OEIridiumCategory_MT:
                    chosen_du = du_bytes
        if not self.args.wanted_du and chosen_du is None:
            chosen_du = du_byte_list[0]
        if chosen_du:
            out_record = OERecord()
            out_record.set_value(SpruceFields().du_single, chosen_du)
            self.success.emit(out_record)
        else:
            self.log.warn("Could not find a suitable DU.")
            self.success.emit(record)


class _SprucePrepBase():
    no_protonate = SpruceParameters().no_protonate
    mutations = SpruceParameters().mutations
    build_sc = SpruceParameters().build_sc
    cap_termini = SpruceParameters().cap_termini
    no_interactions = SpruceParameters().no_interactions
    ligand_name = SpruceParameters().ligand_name
    split_cofactors = SpruceParameters().split_cofactors
    cofactors = SpruceParameters().cofactors
    excipients = SpruceParameters().excipients
    min_atoms = SpruceParameters().min_atoms
    max_atoms = SpruceParameters().max_atoms
    max_residues = SpruceParameters().max_residues
    max_system_atoms = SpruceParameters().max_system_atoms
    site_size = SpruceParameters().site_size

    def begin(self):
        self.prep_opts_args = {'--no-protonate': self.args.no_protonate,
                               '--mutations': self.args.mutations,
                               '--build-sc': self.args.build_sc,
                               '--cap-termini': self.args.cap_termini,
                               '--no-interactions': self.args.no_interactions}

        self.split_opts_args = {'--ligand-name': self.args.ligand_name,
                               '--split-cofactors': self.args.split_cofactors,
                               '--cofactors': self.args.cofactors,
                               '--excipients': self.args.excipients,
                               '--split-cofactors': self.args.split_cofactors,
                               '--min-atoms': self.args.min_atoms,
                               '--max-atoms': self.args.max_atoms,
                               '--max-residues': self.args.max_residues,
                               '--max-system-atoms': self.args.max_system_atoms,
                               '--site-size': self.args.site_size}
        self.split_opts = OESpruceSplitOptions(**self.split_opts_args)
        self.prep_opts = OESprucePrepOptions(**self.prep_opts_args)


class SprucePrepAdvanced(_SprucePrepBase, RecordPortsMixin, ComputeCube):
    title = "Spruce Prep Advanced"
    classification = [["Spruce", "Prep", "Advanced"]]
    tags = ["OpenEye", "snowball", "Spruce"]
    description = "A vector of Spruce-prepped OEDesignUnits is generated from input PDB/MTZ files on the input OERecord."  # noqa

    def process(self, record, port):
        pdb_code = record.get_value(SpruceFields().pdb_code)
        pdb = record.get_value(SpruceFields().pdb_file)
        mol = oechem.OEMol()
        if not ReadProteinFromPDB(pdb.decode(), mol, expand_alts=False):
            self.log.warning("Could not read PDB data for code: {}".format(pdb_code))  # noqa
            self.failure.emit(record)
            return
        mtz_data = None
        if record.has_value(SpruceFields().mtz_file):
            mtz_data = record.get_value(SpruceFields().mtz_file)
        du_ref = None
        if record.has_value(SpruceFields().du_ref_file):
            du_ref_bytes = record.get_value(SpruceFields().du_ref_file)
            du_ref = oechem.OEDesignUnit()
            oechem.OEReadDesignUnitFromBytes(du_ref, du_ref_bytes)
        bu_ref = None
        if record.has_value(SpruceFields().bu_ref_mol):
            bu_ref = record.get_value(SpruceFields().bu_ref_mol)

        du_list, ird_list = create_du_ird(mol, pdb_code, mtz_data, self.split_opts, self.prep_opts, du_ref, bu_ref)  # noqa

        out_record = OERecord()
        out_record.set_value(SpruceFields().du_vec, du_list)
        out_record.add_field(SpruceFields().ird_vec)
        if mtz_data:
            out_record.set_value(SpruceFields().ird_vec, ird_list)

        out_record.add_field(SpruceFields().meta_data)
        if record.has_value(SpruceFields().meta_data):
            out_record.set_value(SpruceFields().meta_data, record.get_value(SpruceFields().meta_data))  # noqa
        out_record.add_field(SpruceFields().mtz_file)
        if record.has_value(SpruceFields().mtz_file):
            out_record.set_value(SpruceFields().mtz_file, record.get_value(SpruceFields().mtz_file))
        out_record.set_value(SpruceFields().pdb_code, record.get_value(SpruceFields().pdb_code))
        out_record.set_value(SpruceFields().pdb_file, record.get_value(SpruceFields().pdb_file))

        self.success.emit(out_record)


class SprucePrepBasic(_SprucePrepBase, RecordPortsMixin, SourceCube):
    title = "Spruce Prep Basic"
    classification = [["Spruce", "Prep", "Basic"]]
    tags = ["OpenEye", "snowball", "Spruce", "Basic"]
    description = "A vector of Spruce-prepped OEDesignUnits is generated from PDB code or set of input PDB/MTZ files."  # noqa

    pdb_code = SpruceParameters().pdb_code
    design_ref = SpruceParameters().design_ref
    biounit_ref = SpruceParameters().biounit_ref

    pdb_file = SpruceParameters().pdb_file
    mtz_file = SpruceParameters().mtz_file

    success = RecordOutputPort('success')
    failure = RecordOutputPort('failure')

    def __iter__(self):

        record = OERecord()

        if self.args.design_ref:
            du_ref = get_du_ref_file(self.args.design_ref)
            if du_ref is None:
                msg = "Could not read the DU ref file: {}. Aborting.".format(self.args.design_ref)
                self.log.warning(msg)  # noqa
                self.failure.emit(record)
                return
        else:
            du_ref = None
        if self.args.biounit_ref:
            bu_ref = get_bu_ref_mol(self.args.biounit_ref)
            if bu_ref is not None:
                msg = " Could not download the BU ref: {}".format(self.args.biounit_ref)
                self.log.warning(msg)  # noqa
                self.failure.emit(record)
                return
        else:
            bu_ref = None

        if self.args.pdb_code:
            if self.args.pdb_file:
                pdb = self.args.pdb_file
            else:
                pdb_url = get_pdb_url(self.args.pdb_code)
                try:
                    pdb = download_pdb(pdb_url)
                except IOError:
                    self.log.info("PDB ({}) Could not download the PDB. Aborting.".format(self.args.pdb_code))  # noqa
                    self.failure.emit(record)
                    return
            pdb_code = self.args.pdb_code
        else:
            if self.args.pdb_file:
                if not os.path.isfile(str(self.args.pdb_file)):
                    if not isinstance(self.args.pdb_file, dict):
                        file = self.args.pdb_file
                    else:
                        file = self.args.pdb_file["file"]
                    resource = APISession.get_resource(File, file)
                    self.temp = TemporaryPath(suffix=resource.name)
                    resource.download_to_file(self.temp.path)
                    self.args.pdb_file = self.temp.path
                if isinstance(self.temp.path, str):
                    with open(self.temp.path, "rb") as pdb_file:
                        pdb = pdb_file.read()
                else:
                    pdb = self.args.pdb_file
            else:
                self.log.info("A PDB file or PDB code is required. Aborting.".format(self.args.pdb_code))  # noqa
                self.failure.emit(record)
                return
        mol = oechem.OEMol()
        if not ReadProteinFromPDB(pdb.decode(), mol, expand_alts=False):
            self.log.warning("Could not read PDB data for code: {}".format(pdb_code))  # noqa
            self.failure.emit(record)
            return
        if not self.args.pdb_code:
            pdb_code = mol.GetTitle() if len(mol.GetTitle()) >= 4 else mol.GetTitle()

        if self.args.mtz_file:
            if not os.path.isfile(str(self.args.mtz_file)):
                if not isinstance(self.args.mtz_file, dict):
                    file = self.args.mtz_file
                else:
                    file = self.args.mtz_file["file"]
                resource = APISession.get_resource(File, file)
                self.temp = TemporaryPath(suffix=resource.name)
                resource.download_to_file(self.temp.path)
                self.args.mtz_file = self.temp.path
            if isinstance(self.temp.path, str):
                with open(self.temp.path, "rb") as mtz_file:
                    mtz = mtz_file.read()
            else:
                mtz = self.args.mtz_file
        else:
            mtz_url = get_mtz_url(pdb_code, self.log.info)
            mtz = None
            if mtz_url:
                mtz = download_mtz(mtz_url, pdb_code, self.log.info)
        du_list, ird_list = create_du_ird(mol, pdb_code, mtz,
                                          self.split_opts, self.prep_opts,
                                          du_ref, bu_ref)

        record.set_value(SpruceFields().du_vec, du_list)
        record.add_field(SpruceFields().ird_vec)
        record.add_field(SpruceFields().mtz_file)
        if mtz:
            record.set_value(SpruceFields().mtz_file, mtz)
            record.set_value(SpruceFields().ird_vec, ird_list)
        record.set_value(SpruceFields().pdb_code, pdb_code)
        record.set_value(SpruceFields().pdb_file, pdb)

        yield record


class ParallelSprucePrepAdvanced(ParallelMixin, SprucePrepAdvanced):
    title = "Parallel " + SprucePrepAdvanced.title
    description = "(Parallel) " + SprucePrepAdvanced.description


class ParallelUrlToFile(ParallelMixin, UrlToFile):
    title = "Parallel " + UrlToFile.title
    description = "(Parallel) " + UrlToFile.description


class ParallelDUtoReceptorDataset(ParallelMixin, DUtoReceptorDataset):
    title = "Parallel " + DUtoReceptorDataset.title
    description = "(Parallel) " + DUtoReceptorDataset.description


class ParallelDUVectorToSingleDU(ParallelMixin, DUVectorToSingleDU):
    title = "Parallel " + DUVectorToSingleDU.title
    description = "(Parallel) " + DUVectorToSingleDU.description


class ParallelDUtoMDDataset(ParallelMixin, DUtoMDDataset):
    title = "Parallel " + DUtoMDDataset.title
    description = "(Parallel) " + DUtoMDDataset.description