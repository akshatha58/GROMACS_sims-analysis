# process_trajs.sh
# ./process_trajs.sh <simulation directory> <simulation time> <GROMACS path>


dir=$1
simtime=$2
GROMACS=$3

cd $dir
mkdir $dir/results_$simtime

trajfile=$dir/md.xtc
trajfileout=$dir/md_centered.xtc

trajfileoutA=$dir/md_chainA.gro
trajfileoutB=$dir/md_chainB.gro
pdbfileoutA=$dir/md_chainA.pdb
pdbfileoutB=$dir/md_chainB.pdb

structfile=$dir/md.tpr
indexfileA=$dir/indexA.ndx
indexfileB=$dir/indexB.ndx
resultsdir=$dir/results_$simtime
energyfile=$dir/md.edr

mkdir $resultsdir/clusters

# Analysis of simulation trajectories
printf "1 0" |  $GROMACS trjconv -f $trajfile -s $structfile -o $trajfileout -pbc whole -center -tu ns 
$GROMACS make_ndx -f $structfile -o $indexfile

# Crosscheck atom ranges for chain A and chain B to insert to index file: recheck
# Chain A: a 1-312
# Chain B: a 313-772 
# name 17 chainA
# name 18 chainB
# q

# Subsequent plotting and analysis is for all trajectories combined

# Get .gro and .pdb files for separate chains according to the index files
printf "17" | $GROMACS trjconv -f $trajfileout -s $structfile -o $trajfileoutA -n $indexfile
printf "18" | $GROMACS trjconv -f $trajfileout -s $structfile -o $trajfileoutB -n $indexfile
printf "17" | $GROMACS trjconv -f $trajfileout -s $structfile -o $pdbfileoutA -n $indexfile
printf "18" | $GROMACS trjconv -f $trajfileout -s $structfile -o $pdbfileoutB -n $indexfile

# # # $GROMACS trjconv -f $trajfileout -s $structfile -o $trajfileoutA -n $indexfile
# # # $GROMACS trjconv -f $trajfileout -s $structfile -o $trajfileoutB -n $indexfile
# # # $GROMACS trjconv -f $trajfileout -s $structfile -o $pdbfileoutA -n $indexfile
# # # $GROMACS trjconv -f $trajfileout -s $structfile -o $pdbfileoutB -n $indexfile

# # # # Crosscheck index file indices before running
printf "17 17" | $GROMACS rms -f $trajfileout -s $structfile -o $resultsdir/rmsd_chainA.xvg -n $indexfile
printf "18 18" | $GROMACS rms -f $trajfileout -s $structfile -o $resultsdir/rmsd_chainB.xvg -n $indexfile
printf "17" | $GROMACS gyrate -f $trajfileout -s $structfile -o $resultsdir/rg_chainA.xvg -n $indexfile
printf "18" | $GROMACS gyrate -f $trajfileout -s $structfile -o $resultsdir/rg_chainB.xvg -n $indexfile
printf "17 17" | $GROMACS rmsf -f $trajfileout -s $structfile -res -o $resultsdir/rmsf_chainA.xvg -n $indexfile
printf "18 18" | $GROMACS rmsf -f $trajfileout -s $structfile -res -o $resultsdir/rmsf_chainB.xvg -n $indexfile

# $GROMACS rms -f $trajfileout -s $structfile -o $resultsdir/rmsd_chainA.xvg -n $indexfile
# $GROMACS rms -f $trajfileout -s $structfile -o $resultsdir/rmsd_chainB.xvg -n $indexfile
# $GROMACS gyrate -f $trajfileout -s $structfile -o $resultsdir/rg_chainA.xvg -n $indexfile
# $GROMACS gyrate -f $trajfileout -s $structfile -o $resultsdir/rg_chainB.xvg -n $indexfile
# $GROMACS rmsf -f $trajfileout -s $structfile -res -o $resultsdir/rmsf_chainA.xvg -n $indexfile
# $GROMACS rmsf -f $trajfileout -s $structfile -res -o $resultsdir/rmsf_chainB.xvg -n $indexfile

# # # # Crosscheck indices before running
printf "11 0" | $GROMACS energy -f $energyfile -o $resultsdir/potential.xvg >> $resultsdir/metrics.txt
printf "12 0" | $GROMACS energy -f $energyfile -o $resultsdir/kinetic.xvg >> $resultsdir/metrics.txt
printf "13 0" | $GROMACS energy -f $energyfile -o $resultsdir/total_energy.xvg >> $resultsdir/metrics.txt
printf "15 0" | $GROMACS energy -f $energyfile -o $resultsdir/temperature.xvg >> $resultsdir/metrics.txt
printf "17 0" | $GROMACS energy -f $energyfile -o $resultsdir/pressure.xvg >> $resultsdir/metrics.txt
printf "23 0" | $GROMACS energy -f $energyfile -o $resultsdir/density.xvg >> $resultsdir/metrics.txt

# # Get secondary structure information with time (single precision calculations only, regardless of GROMACS precision configuration)
$GROMACS dssp -f $dir/md_chainA.pdb -s $dir/md_chainA.gro -o $resultsdir/ss_A.dat -num $resultsdir/ss_data_A.xvg 
$GROMACS dssp -f $dir/md_chainB.pdb -s $dir/md_chainB.gro -o $resultsdir/ss_B.dat -num $resultsdir/ss_data_B.xvg

mkdir $resultsdir/clusters
# printf "17 1" | $GROMACS cluster -f $trajfile -s $structfile -method gromos -cutoff 0.15 -g $resultsdir/clusters/clusters_A.log -cl $resultsdir/clusters/clusters_A.pdb -n $indexfile
# printf "18 1" | $GROMACS cluster -f $trajfile -s $structfile -method gromos -cutoff 0.15 -g $resultsdir/clusters/clusters_B.log -cl $resultsdir/clusters/clusters_B.pdb -n $indexfile
printf "1 1" | $GROMACS cluster -f $trajfileout -s $structfile -method gromos -cutoff 0.15 -g $resultsdir/clusters/clusters.log -cl $resultsdir/clusters/clusters.pdb -n $indexfile
mv $dir/rmsd-clust.xpm $resultsdir/clusters
mv $dir/rmsd-dist.xvg $resultsdir/clusters


$GROMACS make_ndx -f $trajfileoutB -o $indexfileB
$GROMACS make_ndx -f $trajfileoutA -o $indexfileA

# a 1-122: N terminal of chain B, 
# a 349-460: C terminal of chain B, 
# a 1-69: Rf_terminal of chain B
# q

# a 1-161: N terminal of chain A
# q

printf "12 12" | $GROMACS rms -f $pdbfileoutB -s $pdbfileoutB -o $resultsdir/rmsd_Rf_terminal.xvg -n $indexfileB
printf "10 10" | $GROMACS rms -f $pdbfileoutB -s $pdbfileoutB -o $resultsdir/rmsd_ChainB_N.xvg -n $indexfileB
printf "11 11" | $GROMACS rms -f $pdbfileoutB -s $pdbfileoutB -o $resultsdir/rmsd_ChainB_C.xvg -n $indexfileB
printf "10 10" | $GROMACS rms -f $pdbfileoutA -s $pdbfileoutA -o $resultsdir/rmsd_ChainA_N.xvg -n $indexfileA

$GROMACS make_ndx -f $trajfileoutB -o $indexfileB
# $GROMACS make_ndx -f $trajfileoutA -n $indexfileA

# a 5 | a 440 : N and C terminal distances
# a 72 | a 188 : B5 and B13 C-alpha carbon distances (hydrogen bond formation)
# a 391 | a 172 : CT open/close conformations
# a 434 | a 118 : CT open/close conformations

$GROMACS distance -f $pdbfileoutB -s $pdbfileoutB -n $indexfileB -oav $resultsdir/NC_distance.xvg 
$GROMACS distance -f $pdbfileoutB -s $pdbfileoutB -n $indexfileB -oav $resultsdir/B5B13_distance.xvg
$GROMACS distance -f $pdbfileoutB -s $pdbfileoutB -n $indexfileB -oav $resultsdir/B26B12_distance.xvg 
$GROMACS distance -f $pdbfileoutB -s $pdbfileoutB -n $indexfileB -oav $resultsdir/B28B8_distance.xvg

