# truncate_trajs.sh
# 16 February 2025

# Truncates trajectories to 1 ns. Gets .gro, .trr, .xtc files for longer trajectories for analysis
dir=$1
simtime=$2
GROMACS=$3

# mkdir ../$simulation_system/results_$simtime
cd $dir

trajfile=$dir/md.xtc
trajfileout=$dir/../$simtime'_ns'/md.xtc
trrfile=$dir/md.trr
trrfileout=$dir/../$simtime'_ns'/md.trr
structfile=$dir/md.tpr

printf "1 0" |  $GROMACS trjconv -f $trajfile -s $structfile -o $trajfileout -pbc whole -center -tu ns -e $simtime
printf "1 0" |  $GROMACS trjconv -f $trrfile -s $structfile -o $trrfileout -pbc whole -center -tu ns -e $simtime
$GROMACS eneconv -f $dir/md.edr -o $dir/../$simtime'_ns'/md.edr -e 1000

cp $dir/md.tpr $dir/../$simtime'_ns'/md.tpr
cp $dir/md.gro $dir/../$simtime'_ns'/md.gro

