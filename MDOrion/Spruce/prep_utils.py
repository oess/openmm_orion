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
import requests
import tempfile

from openeye import (oechem,
                     oedocking,
                     oespicoli)

from orionclient.session import APISession
from orionclient.types import File
from datarecord.utils import TemporaryPath

from spruce.prep.design_unit_utils import OEIsCovalent
from spruce.prep.du_generator import OEMakeDesignUnits
from spruce.prep.sequence import (OEGetProteinInfoFromASU,
                                  OEMakeMolsFromSequences)
from spruce.prep.iridium_score import (OECalculateDPI,
                                       OECalcIridiumData)
from spruce.utils.surface import OECreateBindingSiteSurface


def get_orion_job_tag():
    return "Job {}".format(os.environ["ORION_JOB_ID"]) if "ORION_JOB_ID" in os.environ else ""  # noqa


def get_pdb_url(code):
    return "http://www.rcsb.org/pdb/files/{}.pdb".format(code)


def get_mtz_url(code, log, log_error=True):
    status_url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/electron_density_statistics/{}".format(code.lower())  # noqa
    r = requests.get(status_url)
    if r.status_code == 404:
        if log_error:
            log('PDB ({}): skipping mtz, no density file found'.format(code))
        return None
    return "http://www.ebi.ac.uk/pdbe/coordinates/files/{}_map.mtz".format(code.lower())  # noqa


def download_pdb(url):
    r = requests.get(url)
    if r.status_code != 200 or r.text.find("is not available") >= 0 or r.text.find("Error") > 0:  # noqa
        raise IOError
    return r.content


def download_mtz(url, code, log):
    if url is None:
        return None
    r = requests.get(url)
    if r.status_code == 200:
        return r.content
    else:
        log('{}: skipping mtz, could not download density file'.format(code))
        return None


def write_mtz(mtz_data, mtz_file):
    with open(mtz_file, 'wb') as ofp:
        ofp.write(mtz_data)


def write_pdb(pdb_data, pdb_file):
    with open(pdb_file, 'w') as ofp:
        ofp.write(pdb_data)


def write_mol_to_bytes(mol):
    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_OEB)
    ofs.openstring()
    oechem.OEWriteMolecule(ofs, mol)
    ofs.close()
    return ofs.GetString()


def read_mol_from_bytes(mol_bytes, use_graphmol=True):
    ifs = oechem.oemolistream()
    ifs.SetFormat(oechem.OEFormat_OEB)
    if not ifs.openstring(mol_bytes):
        return None
    if use_graphmol:
        mol = oechem.OEGraphMol()
    else:
        mol = oechem.OEMol()
    if not oechem.OEReadMolecule(ifs, mol):
        return None
    return mol


def ird_to_str(ird):
    res = "Iridium category: {} LaD: {:.2f} ASaD: {:.2f} DPI: {:.2f} ".format(oechem.OEGetIridiumCategoryName(ird.GetCategory()),
                                                                              ird.GetLigandDensity(),
                                                                              ird.GetAsDensity(),
                                                                              ird.GetDPI())
    res += "POL: {:.1} POAS: {:.1} AltConfs: {:.1} PackRes: {:.1} Excp: {:.1} IrrRFree: {:.1} PossCov: {:.1}".format(     # noqa
            str(ird.GetPartOccupancyLigand()), str(ird.GetPartOccupancyProtein()), str(ird.GetAltConfs()),                     # noqa
            str(ird.GetPackingResidues()), str(ird.GetExcipients()), str(ird.GetIrrationalRfree()), str(ird.GetPossibleCovalent()))  # noqa
    return res


def create_du_ird(mol, pdb_code, mtz, split_opts, prep_opts, du_ref, bu_ref):
    info = OEGetProteinInfoFromASU(mol)
    seq_mols = OEMakeMolsFromSequences(info.actual_sequences)
    du_list = []
    ird_list = []
    for du in OEMakeDesignUnits(mol, pdb_code,
                                split_opts, prep_opts,
                                design_ref=du_ref, biounit_ref=bu_ref,
                                seq_mols=seq_mols):
            if mtz:
                dpi_data = OECalculateDPI(mol)
                with tempfile.TemporaryDirectory() as tmp_dir:
                    mtz_file = os.path.join(tmp_dir, "tmp.mtz")
                    write_mtz(mtz, mtz_file)
                    ird = OECalcIridiumData(du, mtz_file)
                sq = oechem.OEStructureQuality()
                sq.SetIridiumData(ird)
                du.GetImpl().SetStructureQuality(sq)
                ird_list.append(ird_to_str(ird))
            du_list.append(oechem.OEWriteDesignUnitToBytes(du))
    return du_list, ird_list


def _clean_mol_tags(mol):
    tags = ["_is_design_unit_", "__design_unit_title__", "_site_residues_", "_gap_residues_"]  # noqa
    oechem.OEClearStyle(mol)
    oechem.OEClearPDBData(mol)
    oechem.OEDeleteInteractionsHintSerializationData(mol)
    oechem.OEDeleteInteractionsHintSerializationIds(mol)
    for tag in tags:
        if mol.HasData(tag):
            mol.DeleteData(tag)


def _clean_receptor_tags(rec):
    if not oedocking.OEIsReceptor(rec):
        return

    extras = []
    for extra in oedocking.OEReceptorGetExtraMols(rec):
        _clean_mol_tags(extra)
        extras.append(oechem.OEGraphMol(extra))

    oedocking.OEReceptorClearExtraMols(rec)

    for extra in extras:
        oedocking.OEReceptorAddExtraMol(rec, extra)


def DesignUnitToReceptor(du, cache_grids=False, err=print, warn=print):
    receptor = oechem.OEGraphMol()
    mol = oechem.OEMol()
    du.GetComponents(mol, oechem.OEDesignUnitComponents_All ^ oechem.OEDesignUnitComponents_PackingResidues ^ oechem.OEDesignUnitComponents_Ligand)  # noqa
    mol.SetTitle(du.GetTitle())
    ligand = None
    if not du.HasLigand():
        if du.HasSiteResidues():
            surf = OECreateBindingSiteSurface(du, resolution=1.0)
            center = oechem.OEFloatArray(3)
            extents = oechem.OEFloatArray(3)
            oespicoli.OEGetSurfaceCenterAndExtents(surf, center, extents)
            oedocking.OEMakeReceptor(receptor, mol,
                                     center[0], center[1], center[2])
        else:
            err("DesignUnit '{}' does not have a ligand or a site definition".format(du.GetTitle()))  # noqa
            err("Please, re-prepare the DU with a reference that has site or a ligand present")  # noqa
            return None
    else:
        ligand = oechem.OEMol()
        du.GetLigand(ligand)
        if OEIsCovalent(ligand):
            # Delete dummy atoms
            for atom in ligand.GetAtoms(oechem.OEHasAtomicNum(0)):
                ligand.DeleteAtom(atom)
            for atom in mol.GetAtoms(oechem.OEHasAtomicNum(0)):
                mol.DeleteAtom(atom)
        oedocking.OEMakeReceptor(receptor, mol, ligand)
    receptor.SetTitle(du.GetTitle())

    if cache_grids:
        if not du.HasLigand() or not ligand or not ligand.IsValid():
            warn('No bound ligand in {}, skipping scoring grid caching.'.format(du.GetTitle()))  # noqa
        else:
            dock = oedocking.OEDock()
            dock.Initialize(receptor)
            pose = oechem.OEMol()
            dock.DockMultiConformerMolecule(pose, ligand)
            dock.CacheScoringSetup(receptor)

    _clean_receptor_tags(receptor)
    return receptor


def ReadProteinFromPDB(pdb_data, mol, from_stream=True, expand_alts=True):
    ifs = oechem.oemolistream()

    ifs.SetFormat(oechem.OEFormat_PDB)
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default |
                  oechem.OEIFlavor_PDB_DATA |
                  oechem.OEIFlavor_PDB_ALTLOC)

    if from_stream:
        if not ifs.openstring(pdb_data):
            oechem.OEThrow.Fatal("Unable to open the pdb data for reading.")
    else:
        if not ifs.open(pdb_data):
            oechem.OEThrow.Fatal("Unable to open PDB data for reading.")

    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from PDB data.")
    ifs.close()

    if expand_alts:
        fact = oechem.OEAltLocationFactory(mol)
        mol.Clear()
        fact.MakePrimaryAltMol(mol)
    return bool(mol)


def read_du(file_path):
    du = oechem.OEDesignUnit()
    if oechem.OEReadDesignUnit(file_path, du):
        return du
    else:
        return None


def get_bu_ref_mol(bu_path):
    if not isinstance(bu_path, dict):
        file_path = bu_path
    else:
        file_path = bu_path["file"]
    if not os.path.isfile(file_path):
        resource = APISession.get_resource(File, file_path)
        temp = TemporaryPath(suffix=resource.name)
        resource.download_to_file(temp.path)
        file_path = temp.path
        bu = read_du(file_path)
        if not bu:
            mol = oechem.OEMol()
            if ReadProteinFromPDB(file_path, mol, from_stream=False):
                return mol
            else:
                return None
    else:
        mol = oechem.OEMol()
        if ReadProteinFromPDB(file_path, mol, from_stream=False):
            return mol
        else:
            return None


def get_du_ref_file(du_path):
    if not isinstance(du_path, dict):
        file_path = du_path
    else:
        file_path = du_path["file"]
    if not os.path.isfile(file_path):
        resource = APISession.get_resource(File, file_path)
        temp = TemporaryPath(suffix=resource.name)
        resource.download_to_file(temp.path)
        file_path = temp.path
    du = read_du(file_path)
    if not du:
        mol = oechem.OEMol()
        if ReadProteinFromPDB(file_path, mol, from_stream=False):
            du = oechem.OEDesignUnit()
            oechem.OEUpdateDesignUnit(du, mol, oechem.OEDesignUnitComponents_Protein)
            return du if du.IsValid() else None
        else:
            return None
    else:
        return du
