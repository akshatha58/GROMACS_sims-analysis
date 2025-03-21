# automate_sims.sh

# ./automate_sims.sh <GROMACS path> <simulation directory> <precision> <PDB ID>

GROMACS=$1
simulation_system=$2
precision=$3
pdb=$4

# Running simulations 

if [ "$3" == "Double" ]
then
    chmod +x Run_Gromacs_doubleprecision.sh
    for i in {0..16}
    do
    ./Run_Gromacs_doubleprecision.sh $GROMACS $simulation_system $i $pdb
    done

elif [ "$3" == "Single" ]
then 
    chmod +x Run_Gromacs_singleprecision.sh
    for i in {0..16}
    do
    ./Run_Gromacs_singleprecision.sh $GROMACS $simulation_system $i $pdb
    done     
fi 

