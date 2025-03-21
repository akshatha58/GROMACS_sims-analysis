# automate_sims.sh
# 11 February 2025

# Directory names:
# gromacs_dp, gromacs_sp, calligo_dp, calligo_sp (all simulations are for merged chains from now on)

GROMACS=$1
simulation_system=$2
precision=$3
pdb=$4


# Running simulations 

if [ "$3" == "Double" ]
then
    echo "Double"
    chmod +x Run_Gromacs_doubleprecision.sh
    for i in {1..16}
    do
    ./Run_Gromacs_doubleprecision.sh $GROMACS $simulation_system $i $pdb
    done

elif [ "$3" == "Single" ]
then 
    echo "Single"
    chmod +x Run_Gromacs_singleprecision.sh
    for i in {1..16}
    do
    ./Run_Gromacs_singleprecision.sh $GROMACS $simulation_system $i $pdb
    done     
fi 

