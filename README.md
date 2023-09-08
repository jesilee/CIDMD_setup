# CIDMD_setup
#### Authors: Jesi Lee, Lee-Ping Wang


## General Description
* These scripts set up CIDMD (collision-indueced dissociation molecular dynamics) simulations in TeraChem.
* Given a molecular ion (protonated molecule), CIDMD predicts a theoretical CID-MS/MS spectrum by automatically placing a collider gas argon atom and accelerating the collider towards the molecular ion using NEV _ab initio_ MD simulations.
* To run the scripts: ```bash CIDMD_setup.com```


## Organization
For the scripts to work properly, the four files, `CIDMD_setup.com`, `CIDMD_setup.py`, `CIDMD_preprun.py`, and `lebedev.py`, should be located in the current working directory.
`CIDMD_setup.com` contains all variables that users can set for `CIDMD_setup.py` to take, executes the rest of scripts, and creates number (XXX) of input files in the subfolder structures: cid/calcs/XXX/chunk_0000/ (where XXX ranges from 000 to the user-defined three digit number). 
To see all variables, see the *Inputs* section below.
Each resulting subfolder contains two files: `start.xyz` and `vel.xyz`.
* `start.xyz` contains the coordinates for the molecular ion user-provided and the collider gas argon atom cooridiates located an user-defined distance away from the molecular ion.
* `vel.xyz` contains the velocity for the molecular ion user-provided (default=0) and the user-defined velocity (in AKMA unit) of the collider gas argon atom.


## Inputs
* an optimized molecular ion structure in xyz file format.
* variables for `CIDMD_setup.py` (please see below for the variable options from the help section: python3 CIDMD_setup.py -h)
```
usage: CIDMD_setup.py [-h] [--extent EXTENT] [--cutoff CUTOFF]
                            [--set_distance SET_DISTANCE]
                            [--n_points N_POINTS] [--seed SEED]
                            [--deterministic] [--cartgrid CARTGRID]
                            [--lebedev LEBEDEV] [--v_scale V_SCALE] [-v]
                            xyz_in

positional arguments:
  xyz_in                REQUIRED positional argument: input .xyz structure of
                        molecular ion

optional arguments:
  -h, --help            show this help message and exit
  --extent EXTENT       Half the edge-length of cubic box in Angstrom in which
                        initial positions are generated.Should be larger than
                        the spatial extent of the molecule. (default: 50.0)
  --cutoff CUTOFF       Keep only the collider initial positions whose paths
                        come to within this distance of any atom on the
                        molecular ion (in Angstrom) (default: 2.0)
  --set_distance SET_DISTANCE
                        Atomic radius of molecular envelope surface on which
                        the collider will be placed (in Angstrom). (default:
                        5.0)
  --n_points N_POINTS   Number of initial conditions desired when using
                        randomized algorithm. (default: 50)
  --seed SEED           Number of initial conditions desired when using
                        randomized algorithm. (default: 2021)
  --deterministic       If set, generate initial positions using deterministic
                        algorithm instead of default randomized algorithm.
                        (default: False)
  --cartgrid CARTGRID   Number of Cartesian grid points for deterministic
                        algorithm. Default is set to generate 51 points as
                        -50, -48, -46, ..., 48, 50. (default: 51)
  --lebedev LEBEDEV     Lebedev grid order for deterministic algorithm (use
                        odd number in range [3,31]; higher orders not
                        recommended) (default: 5)
  --v_scale V_SCALE
  -v, --verbose         Turn on additional output for debugging. (default:
                        False)

```  
## Outputs
* a user-defined number of CIDMD input files, `start.xyz` and `vel.xyz`, varing collisional geometry
* a user-defined number of collision trajectories that can be studied for fragmentation reaction pathways
* further analyses can be done using this repositories:
  `github.com/jesilee/CIDMD_analysis` performs analysis of CIDMD results, and
  `github.com/jesilee/CIDMD_compare` compares the quality of CIDMD prediction against the reference mass spectra provided.

## Dependencies
* python 3.6
* numpy 

* these files from this repo should be in the current working directory: 
 `CIDMD_setup.com`
 `CIDMD_setup.py`
 `CIDMD_preprun.py`
 `lebedev.py`


## Important notes
To start CIDMD, Folder organization is very important. Here is the instruction to set up CIDMD with this repository:

1) start by creating a folder `Molecule` where molecule is the name of user-defined molecule.
   example: Psilcin
   
2) create a text file name `mol_info.in` that contains information about the molecule.
   example:
```
###__ mol_info __###
mol_name = Psiocin
mol_id   = PS
mol_mf   = C12H16N2O
mol_mw   = 204
mol_xw   = 204.268
```

3) git clone this current repository:
  `git clone https://github.com/jesilee/CIDMD_setup`

4) check these files are in the current working directory: 
  `CIDMD_setup.com`, `CIDMD_setup.py`, `CIDMD_preprun.py`, and `lebedev.py`

5) decid the variable options for CIDMD (please see *Inputs* section) and change the variables in `CIDMD_setup.com`
   
6) execute `bash CIDMD_setup.com`.
   
7) run these trajectories with TeraChem
   
8) process CIDMD trajectories using Learnreaction.py
    
9) analyze CIDMD trajectories using these repositories:
  `git clone https://github.com/jesilee/CIDMD_analysis`
  `git clone https://github.com/jesilee/CIDMD_compare`  




