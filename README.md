# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `LigPrepCubes/` - Cubes for preparing molecules
  * `OEOmegaConfGen` - Use OE's omega cubes to generate MCMols
  * `FREDDocking` - Dock MCMols using FRED to a prepared receptor
  * `SMIRFFParameterization` - Parametrizes molecule with SMIRFF
  * `SetIDTagfromTitle` - Attaches an idtag from the molecule's Title or random ID.
  * `OEBSinkCube` - Custom cube to write out compressed MCMols
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - set up protein:ligand complex and emit OpenMM System
  * `OpenMMSimulation` - run OpenMM simulation.

## Workfloes

* `floes/openmm_complex_setup.py` -  setup the protein:ligand complex.
* `floes/openmm_setup_md.py` - setup the protein:ligand complex and run 1000 MD steps.
* `floes/openmm_md.py` - run MD simulation from a prepared complex.
* `floes/openmm_continue.py` - restart MD simulation.
* `floes/ligprep_floes.py` - parse smiles, dock, and parameterize molecules.
* `floes/setupmd_floe.py` - parse smiles, prepare protein:ligand complex, run MD simulation.

## Local Installation
### MacOS 10.12 Sierra
```bash
git clone git@github.com:openeye-private/openmm_orion.git
cd openmm_orion

#Create a new local conda environment and install dependencies
conda create -c omnia -c omnia/label/dev -n openmm_mac python=3.5 openmm==7.0.1 openmoltools==0.7.4 ambermini==16.16.0 smarty==0.1.4 parmed==2.7.1
source activate openmm_mac

#Install PDBFixer (pinned to specific commit for now)
pip install -e git+https://github.com/pandegroup/pdbfixer.git@5ed0d2550b156961ae4de900f33ae6c6120faea7#egg=pdbfixer

#Install older stable version of OE Toolkits
pip install OpenEye-toolkits-python3-osx-10.11-x64-2016.10.1.tar.gz

#Modify the OpenEye-floe installation requirements
tar -xvzf OpenEye-floe-0.2.127.tar.gz
cd OpenEye-floe-0.2.115
#Change line to install_requires=['requests']
vi setup.py
#Install the OpenEye-floe package
python setup.py develop

#Back into the main folder, Install the main OpenMM Orion Floes
cd ..
python setup.py develop
#Run the Example below
```

### Linux
```bash
git clone git@github.com:openeye-private/openmm_orion.git
cd openmm_orion

#Create a new local conda environment and install dependencies
conda create -c omnia -c omnia/label/dev -n openmm_linux python=3.5 openmm==7.0.1 openmoltools==0.7.4 ambermini==16.16.0 smarty==0.1.4 parmed==2.7.1 pdbfixer-dev
source activate openmm_linux

#Install the OpenEye-floe package
pip install OpenEye-floe-0.2.127.tar.gz

#Install the main OpenMM Orion Floes
python setup.py develop
#Run the Example below
```

## Example
### Starting from SMILES
```bash
# Starting from SMILES, generate conformers with OMEGA
# Dock using FRED, and parameterize with SMIRFF parameters.
python floes/ligprep_floe.py --ifs input/test_smiles.ism --receptor input/test-receptor.oeb.gz --ffxml input/smirff99Frosst.ffxml

# Does everything above and prepares the complex, run 1000 MD steps.
python floes/setupmd_floe.py --ifs input/test_smiles.ism --receptor input/test-receptor.oeb.gz --ffxml input/smirff99Frosst.ffxml --protein input/receptor-fixed.pdb
```

### Starting from a PDB
```bash
# Setup protein-ligand complex
python floes/openmm_complex_setup.py --ligand input/toluene.pdb --protein input/T4-protein.pdb --ffxml input/smirff99Frosst.ffxml

# Setup protein-ligand complex and Run short MD simulation
python floes/openmm_setup_md.py --ligand input/toluene.pdb --protein input/T4-protein.pdb --ffxml input/smirff99Frosst.ffxml

# Above with BLUES
python floes/blues_floe.py --ligand input/toluene.pdb --protein input/T4-protein.pdb --ffxml input/smirff99Frosst.ffxml --steps 5000

# Run a MD Simulation given a complex oeb.gz
python floes/openmm_md.py --ifs input/9PC1X-complex.oeb.gz

# Restart a simulation from saved state
python floes/openmm_continue.py --ifs output/9PC1X-simulation.oeb.gz
```
