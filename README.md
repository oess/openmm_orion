# OpenMM cubes and Floes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
  * `PlatformTestCube`- Checks available OpenMM Platforms
  * `BenchmarkCube` - Benchmark available OpenMM Platforms
* `LigPrepCubes/` - Cubes for preparing molecules
  * `LigChargeCube`- Charges ligands by using the ELF10 charge method
  * `FREDDocking` - Dock MCMols using FRED to a prepared receptor
* `ComplexPrepCubes`
  * `HydrationCube` - Cube for hydrate a molecular system
  * `SolvationCube` - Cube for Solvate a molecular system in a given mixture solvent.
  * `ComplexPrep` - Assembles complex made of solvated system and ligands
  * `ForceFieldPrep` - Parameterizes the system given a forcefield.
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMminimizeCube` - Minimize the protein:ligand complex.
  * `OpenMMnvtCube` - NVT simulation of the protein:ligand complex
  * `OpenMMnptCube`- NPT simulation of the protein:ligand complex.
* `YankCubes/` - YANK cubes
  * `YankHydrationCube` - YANK hydration free energy calculations
  * `YankBindingCube` - YANK absolute binding free energy calculations
  * `YankSolvationFECube` - YANK solvation free energy calculations
  * `YankBindingFECube` - YANK absolute binding free energy calculations
   

## Workfloes
* Testing Floes:
  * `floes/platformTest.py` - Check available OpenMM Platforms
  * `floes/openmm_benchmarking.py` - Performs Benchmarking upon all available Platforms.
* [OpenMM](https://github.com/pandegroup/openmm) Floes:
  * `floes/openmm_complex_prep.py` - Assemble the complex starting from a set of ligands and a given receptor
  * `floes/openmm_complex_prep_min.py` - Assemble and minimize the complex starting from a set of ligands and a given receptor
  * `floes/openmm_MDminimize.py` - Minimize an OpenMM-ready solvated complex
  * `floes/openmm_MDnpt.py` - NPT simulation of an OpenMM-ready solvated complex
  * `floes/openmm_MDnvt.py` - NVT simulation of an OpenMM-ready solvated complex
  * `floes/openmm_MDprep.py` - Set up an OpenMM complex then minimize, warm up and equilibrate a system by using three equilibration stages
  * `floes/openmm_MDprod.py` - Run an unrestrained NPT simulation at 300K and 1atm
  * `floes/openmm_MDprep_prod.py` - Set up an OpenMM complex then minimize, warm up and equilibrate a system by using three equilibration stages. 
  Finally a 2ns production simulation is performed     
* [YANK](https://github.com/choderalab/yank) Floes
  * `floes/yank_hydration.py` - Compute small molecule hydration free energies using YANK.
  * `floes/yank_binding.py` - Compute small molecule absolute binding free energies using YANK.
  * `floes/solvation_free_energy.py` - Compute small molecule solvation free energies using YANK.
  * `floes/binding_free_energy.py` - Compute small molecule absolute binding free energies using YANK.
 
## Local Installation
```bash
git clone https://github.com/oess/openmm_orion.git
cd openmm_orion

#Create a new local conda environment and install dependencies
conda create -n openmmorion python=3.5
source activate openmmorion
conda install -c omnia -c omnia/label/dev -c mobleylab -c OpenEye/label/Orion -c conda-forge openmm==7.1.1 openmoltools==0.8.1 ambermini==16.16.0 parmed==2.7.3 pdbfixer==1.4 openforcefield==0.0.2 smirff99frosst==1.0.5 alchemy==1.2.3 yank==0.18.0 oeommtools pymbar==3.0.3 networkx==1.11 

#Install the OpenEye-floe package and toolkits
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
pip install --pre --extra-index-url https://pypi.anaconda.org/OpenEye/channel/beta/simple OpenEye-oenotebook
pip install OpenEye-floe-0.2.181.tar.gz

#Install the main OpenMM Orion Floes
python setup.py develop

# Run the tests.
py.test -v -s PlatformTestCubes
py.test -v -s LigPrepCubes
py.test -v -s ComplexPrepCubes
py.test -v -s -m "not slow" OpenMMCubes
```
