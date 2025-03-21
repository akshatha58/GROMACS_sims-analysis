# GROMACS_sims-analysis
Code to automate GROMACS simulations, trajectory processing and analysis of certain structural and thermodynamic properties.

## Overall directory structure :

```bash
├── GROMACS_sims-analysis
│   ├── Base_Files                             > (.pdb and .mdp input files) 
│   ├── user_inputs                            > (user inputs for analysis_codes)
│   ├── Simulation directories
│   │   ├── Simulation directory <i>
│   │   │   ├── results_<simtime>              > All results in .xvg/.dat format
│   │   │___│___├── clusters                   > .log and .pdb files from clustering analysis
│   │   ├── plots_trial                        > .png images for general plots
│   │   │   ├── distances                      > .png files for distance measures
│   │   │   ├── smooth                         > .png files: smoothened plots
│   │___│___├── terminal_RMSDs                 > .png files: terminal RMSD information
│___├── analysis_codes                         > (all simulation/processing/analysis codes)
```

## Main run code: 

The code sequentially runs and processes simulations for a list of simulation replicates or precisions (single or double). Can also extend it to automate runs for multiple simulation times or multiple starting structures. The analysis and plotting is then done for all simulations in the list together, to make it easier to visualise and compare.


IMPORTANT NOTES: 
1. Change user inputs at the beginning of ./main.sh before running (Lines 4 to 10).
2. Also change GROMACS paths in lines 34 and 39 for single and double precision.
3. For GROMACS user inputs (such as force field choice, disulphide linking, and index file calls), cross check parameter numbers from `user_inputs/parameter_numbers.txt`. For index file calls, crosscheck corresponding atom numbers from the `md.gro` file.
4. If a simulation needs to be stopped and continued from a checkpoint, execute the following command:
   ```
   # GROMACS=gmx_file_path
   # simulation_system=simulation_directory_path
   # pdb=PDB_ID
   
   ./Run_Gromacs_doubleprecision.sh $GROMACS $simulation_system 17 $pdb
   # OR
   ./Run_Gromacs_singleprecision.sh $GROMACS $simulation_system 17 $pdb
   ```
   

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

## Sequence of execution in analysis_codes:

-iterate over each directory-

* Creation of folders according to the standard tree format above
* Transfer of base .pdb and .mdp files to the simulation directory
* Simulation runs (`automate_sims.sh`)
* Has a script for truncating trajectories (to a lower simulation time than runtime) (`truncate_trajs.sh`): uncomment if needed
* Processing of trajectories and calculations for thermodynamic/structural properties (`process_trajs.sh`)
* Processing of cluster.pdb to add chain information for visualising in Pymol (`process_pdb.py`)

-end of iteration-

* Analysis data across all simulations, generation of plots (`summary_plots.py`)
  - Plots and statistics for energies, pressure, density, temperature
  - Plots and statistics for RMSD, RMSF, Rg for each chain
  - Plots for secondary structure counts for each chain of insulin(with time)
  - Plots for RMSDs of terminal regions of each chain of insulin
  - Plots for all important pairwise CA distances that signify conformational transitions in insulin

