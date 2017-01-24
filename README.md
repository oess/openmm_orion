# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - set up protein:ligand complex and emit OpenMM System
  * `OpenMMSimulation` - run OpenMM simulation.

## Workfloes

* `floes/example_floe.py` - report available OpenMM Platforms
* `floes/openmm_setup_md.py` - setup the protein:ligand complex and run 1000 MD steps.


## Local Installation
### MacOS 10.12 Sierra
```bash
#Create a new local conda environment and install dependencies
conda create -c omnia -c omnia/label/dev -n openmm_mac python=3.5 openmm==7.0.1 openmoltools==0.7.4 ambermini==16.16.0 smarty==0.1.4 parmed==2.7.1
source activate openmm_mac

#Install PDBFixer (pinned to specific commit for now)
pip install -e git+https://github.com/pandegroup/pdbfixer.git@5ed0d2550b156961ae4de900f33ae6c6120faea7#egg=pdbfixer

#Install older stable version of OE Toolkits
pip install OpenEye-toolkits-python3-osx-10.11-x64-2016.10.1.tar.gz

#Modify the OpenEye-floe installation requirements
tar -xvzf OpenEye-floe-0.2.115.tar.gz
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
#Create a new local conda environment and install dependencies
conda create -c omnia -c omnia/label/dev -n openmm_linux python=3.5 openmm==7.0.1 openmoltools==0.7.4 ambermini==16.16.0 smarty==0.1.4 parmed==2.7.1 pdbfixer-dev
source activate openmm_linux

#Install the OpenEye-floe package
pip install OpenEye-floe-0.2.115.tar.gz

#Install the main OpenMM Orion Floes
python setup.py develop
#Run the Example below
```

## Example
```bash
# Setup protein-ligand complex and Run short MD simulation
python floes/openmm_setup_md.py --protein OpenMMCubes/tests/input/T4-protein.pdb --ligand OpenMMCubes/tests/input/smirff_mol.oeb.gz

```
