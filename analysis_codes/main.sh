# main.sh

# simulation_systems=('/home/dynamics/akshatha/Calligo/previous_simfiles/1xda_sp_merged_chain' '/home/dynamics/akshatha/Calligo/previous_simfiles/1xda_dp_merged_chain' '/home/dynamics/akshatha/Calligo/Calligo_sp' '/home/dynamics/akshatha/Calligo/Calligo_dp_merged' '/home/dynamics/akshatha/Calligo/gromacs_sp' '/home/dynamics/akshatha/Calligo/gromacs_dp')

# Inputs: change all inputs in this file
maindir=dummy_files
results=plots_trial
simulation_systems=('/home/dynamics/akshatha/Calligo/dummy_files/gromacs_sp_trial' '/home/dynamics/akshatha/Calligo/dummy_files/gromacs_sp_trial2')
precision=('Single' 'Single')
simtime=1
smooth=1
pdb=1xda

mkdir ../$maindir
basedir=$(pwd)
resultsdir=../$maindir/$results


# Codeblock 1: Running simulations, processing trajectories, and obtaining data from simulations
for i in {0..1}; do

    echo "Creating folders..." 
    echo ${simulation_systems[$i]}

    mkdir ${simulation_systems[$i]}
    mkdir ${simulation_systems[$i]}/$simtime"_ns"

    cp $basedir/../Base_Files/* ${simulation_systems[$i]}/$simtime"_ns"
    cd ${simulation_systems[$i]}/$simtime"_ns"

    dirname=$(basename ${simulation_systems[$i]})
    echo "Running " $dirname "in " ${precision[$i]} "precision..."

    if [ ${precision[$i]} == "Double" ]
    then
    GROMACS=/home/engineering/gromacs-2024.4/build/bin/gmx_d # CHANGE PATH

    elif [ ${precision[$i]} == "Single" ]
    then
    echo "single"
    GROMACS=/usr/local/gromacs/bin/gmx # CHANGE PATH
    fi

    cd $basedir
    ./automate_sims.sh $GROMACS ${simulation_systems[$i]}/$simtime"_ns" ${precision[$i]} $pdb
    
    cd $basedir
    echo "Processing trajectories and obtaining structural and thermodynamic data..."
    # ./truncate_trajs.sh ${simulation_systems[$i]}/$simtime"_ns" $simtime $GROMACS
    ./process_trajs.sh ${simulation_systems[$i]}/$simtime"_ns" $simtime $GROMACS
    python3 process_pdb.py ${simulation_systems[$i]}/$simtime"_ns/results_$simtime/clusters/clusters.pdb"

done

# Codeblock 2: Analysing and plotting data obtained from simulation results (across all trajectories listed in simulation_systems)

echo "Analysing all data..."
cd $basedir
mkdir $resultsdir
mkdir $resultsdir/distances
mkdir $resultsdir/smooth
mkdir $resultsdir/terminal_RMSDs
python3 summary_plots.py "$simtime" "$smooth" "$resultsdir" "${simulation_systems[@]}" > $resultsdir/output.txt



