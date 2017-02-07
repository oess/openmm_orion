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

* `floes/openmm_complex-setup.py` - Setup the protein:ligand complex from **PDBs**.
* `floes/openmm_setup-md.py` - Setup the protein:ligand complex and run MD simulation from **PDBs**
* `floes/openmm_md.py` - Run MD simulation from a prepared *complex.oeb.gz*
* `floes/openmm_restart.py` - Restart MD simulation from a saved state in the *simulation.oeb.gz*.
* `floes/openmm_mdblues.py` - Run MD simulation and MD+BLUES from a *complex.oeb.gz*.
* `floes/smiles_ligprep.py` - Parse smiles, dock and parameterize the molecules.
* `floes/smiles_complex-setup.py` - Parse smiles and setup the protein:ligand complex.
* `floes/smiles_setup-md.py` - Parse smiles, prepare complex, and run MD simulation.
* `floes/smiles_mdblues.py` - Parse smiles, prepare+equilibrate complex, and run MD+BLUES simulation.

## Local Installation
```bash
git clone git@github.com:openeye-private/openmm_orion.git
cd openmm_orion

#Create a new local conda environment and install dependencies
conda create -c omnia -c omnia/label/dev -n dev python=3.5 openmm==7.1.0rc1 openmoltools==0.7.4 ambermini==16.16.0 smarty==0.1.4 parmed==2.7.1 pdbfixer==1.4
source activate dev

#Install the OpenEye-floe package
pip install OpenEye-floe-0.2.127.tar.gz

#Install the main OpenMM Orion Floes
python setup.py develop
```
#### MacOS 10.12 Sierra
If you're running on MacOS 10.12 Sierra, you may have to install
older version of OE Toolkits.
```bash
pip install OpenEye-toolkits-python3-osx-10.11-x64-2016.10.1.tar.gz

#Modify the OpenEye-floe installation requirements
tar -xvzf OpenEye-floe-0.2.127.tar.gz && cd OpenEye-floe-0.2.127

#Change line to install_requires=['requests']
vi setup.py

#Install the OpenEye-floe package from within it's directory.
(OpenEye-floe.0.2.127/) $> python setup.py develop

#Back into the main folder, Install the main OpenMM Orion Floes
cd .. && python setup.py develop
```

## Examples
### Starting from PDBs
#### [SetupOpenMMComplex]: Setup the protein:ligand complex.
*Assuming input ligand PDB is in docked position* Setup the protein-ligand complex.
```
python floes/openmm_complex-setup.py --ligand input/toluene.pdb --protein input/T4-protein.pdb --ffxml input/smirff99Frosst.ffxml
```

#### [SetupOpenMMSimulation]: Setup and simulate the protein:ligand complex.
Does the same above and run 10000 MD steps.
```
python floes/openmm_setup-md.py --ligand input/toluene.pdb --protein input/T4-protein.pdb --ffxml input/smirff99Frosst.ffxml --steps 10000
```

### Using the prepared complex.oeb.gz files
#### [RunOpenMMSimulation] Run a MD Simulation given a complex.oeb.gz
If you want to just run MD from the prepared complex:
```
python floes/openmm_md.py --complex input/9PC1X-complex.oeb.gz --steps 10000
```

#### [RestartOpenMMSimulation] Restart MD Simulation from the saved state in the simulation.oeb.gz
To restart the MD simulation from the previous run:
```
python floes/openmm_restart.py --complex output/9PC1X-simulation.oeb.gz --steps 10000
```

#### [RunOpenMMBLUES] Run an equil MD simulation and MD+Blues
To run MD and then MD+BLUES:
```
python floes/openmm_mdblues.py --complex input/9PC1X-complex.oeb.gz --steps 5000 --mdsteps 5000 --ncsteps 2500
```

### Starting from SMILES
#### [SmilesLigPrep]: Generate conformers, dock and parameterize molecules
Starting from a SMILE string, generate conformers with OMEGA, dock to a
prepared receptor using FRED, and parameterize the molecules. Writes out the
SMIRFF molecules in their docked poses.
```
python floes/smiles_ligprep.py --ligand input/test_smiles.ism --receptor input/test-receptor.oeb.gz --ffxml input/smirff99Frosst.ffxml
```

#### [SmilesComplexPrep]: Setup the protein:ligand complexes
Does the same as above and then prepares the complex from a PDB of the receptor.
Writes out the protein:ligand complex.
```
python floes/smiles_complex-setup.py --ligand input/test_smiles.ism --receptor input/test-receptor.oeb.gz --ffxml input/smirff99Frosst.ffxml --protein input/receptor-fixed.pdb
```

#### [SmilesSimulation]: Setup and prepare the MD simulation.
Does all the preparation steps above and runs the MD simulation:
```
python floes/smiles_setup-md.py --ligand input/test_smiles.ism --receptor input/test-receptor.oeb.gz --ffxml input/smirff99Frosst.ffxml --protein input/receptor-fixed.pdb --steps 10000
```

#### [SmilesMDBLUES] Running all the cubes together.
Prepare the systems, equil. 5000 MD steps, and run 5000 MD+2500 BLUES steps.
```
python floes/smiles_mdblues.py --ligand input/test_smiles.ism --receptor input/test-receptor.oeb.gz --ffxml input/smirff99Frosst.ffxml --protein input/receptor-fixed.pdb --steps 5000 --mdsteps 5000 --ncsteps 2500
```
