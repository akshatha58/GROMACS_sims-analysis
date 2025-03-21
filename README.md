# GROMACS_sims-analysis
Code to automate GROMACS simulations, trajectory processing and analysis of certain structural and thermodynamic properties.

## Overall directory structure :

```bash
├── GROMACS_sims-analysis
│   ├── analysis_codes                         > (all simulation/processing/analysis codes)
│   ├── Base_Files                             > (.pdb and .mdp input files) 
│   ├── user_inputs                            > (user inputs for analysis_codes)
│   ├── Simulation directories
│   │   ├── Simulation directory <i>
│   │   │   ├── results_<simtime>              > All results in .xvg/.dat format
│   │   │___│___├── clusters                   > .log and .pdb files from clustering analysis
│   │   ├── plots_trial                        > .png images for general plots
│   │   │   ├── distances                      > .png files for distance measures
│   │   │   ├── smooth                         > .png files: smoothened plots
│___│___│___├── terminal_RMSDs                 > .png files: terminal RMSD information
```

## Main run code: 

The code sequentially runs and processes simulations for a list of simulation replicates or precisions (single or double). Can also extend it to automate runs for multiple simulation times or multiple starting structures. The analysis and plotting is then done for all simulations in the list together, to make it easier to visualise and compare.


NOTE: 
1. Change user inputs at the beginning of ./main.sh before running (Lines 4 to 10).
2. Also change GROMACS paths in lines 34 and 39 for single and double precision. 

```
maindir=trial_files   # "Simulation directories" mentioned in the tree above
results=plots_trial   # "plots_trial": stores all plot information

# Full paths for Simulation directory <i> (i: total number of simulations run)
simulation_systems=('../trial_files/gromacs_sp_trial' '../trial_files/gromacs_sp_trial2')

precision=('Single' 'Single') # Options are 'Single' or 'Double'

simtime=1 # Units: ns. Remember to change the corresponding entry in the md.mdp file if there is a change here.
smooth=1  # 1 if Gaussian smoothening of data is needed. 0 if you need raw data to be plotted
pdb=1xda  # PDB ID (must match the name of the .pdb file used as the starting structure for simulations)
```

Executing the code:
```
./main.sh
```
