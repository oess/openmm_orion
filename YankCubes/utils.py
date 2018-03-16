import os
import numpy as np
from simtk import unit
import yaml
from yank.analyze import get_analyzer
from floe.api.orion import in_orion
from big_storage import LargeFile


def analyze_directory(source_directory):
    """
    This Function has been copied and adapted from the Yank ver 0.17.0 source code
    (yank.analyse.analyze_directory)

    Analyze contents of store files to compute free energy differences.

    This function is needed to preserve the old auto-analysis style of YANK. What it exactly does can be refined
    when more analyzers and simulations are made available. For now this function exposes the API.

    Parameters
    ----------
    source_directory : string
        The location of the simulation storage files.

    """
    analysis_script_path = os.path.join(source_directory, 'analysis.yaml')
    if not os.path.isfile(analysis_script_path):
        err_msg = 'Cannot find analysis.yaml script in {}'.format(source_directory)
        raise RuntimeError(err_msg)
    with open(analysis_script_path, 'r') as f:
        analysis = yaml.load(f)

    data = dict()
    for phase_name, sign in analysis:
        phase_path = os.path.join(source_directory, phase_name + '.nc')
        phase = get_analyzer(phase_path)
        data[phase_name] = phase.analyze_phase()
        kT = phase.kT

    # Compute free energy and enthalpy
    DeltaF = 0.0
    dDeltaF = 0.0
    DeltaH = 0.0
    dDeltaH = 0.0
    for phase_name, sign in analysis:
        DeltaF -= sign * (data[phase_name]['DeltaF'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaF += data[phase_name]['dDeltaF'] ** 2
        DeltaH -= sign * (data[phase_name]['DeltaH'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaH += data[phase_name]['dDeltaH'] ** 2
    dDeltaF = np.sqrt(dDeltaF)
    dDeltaH = np.sqrt(dDeltaH)

    DeltaF = DeltaF * kT / unit.kilocalories_per_mole
    dDeltaF = dDeltaF * kT / unit.kilocalories_per_mole
    DeltaH = DeltaH * kT / unit.kilocalories_per_mole
    dDeltaH = dDeltaH * kT / unit.kilocalories_per_mole

    return DeltaF, dDeltaF, DeltaH, dDeltaH


def upload(filename):

    file_id = filename

    if in_orion():
        file_id = LargeFile(filename)

    return file_id


def download(file_id):

    filename = file_id

    if in_orion():
        filename = file_id.retrieve()
        file_id.delete()

    return filename