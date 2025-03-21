#!/bin/bash

GROMACS=$1
sim_directory=$2
pdb=$4

cd $sim_directory

if [ "$3" == "0" ]
then
grep -v HETATM $pdb.pdb > $pdb"_new".pdb
elif [ "$3" == "1" ]
then
	 $GROMACS pdb2gmx -f $pdb"_new".pdb -o $pdb"_processed".gro -water spce -ignh -ter -ss -merge all
elif [ "$3" == "2" ]
then
	 $GROMACS editconf -f $pdb"_processed".gro -o $pdb"_box".gro -c -d 2.0 -bt dodecahedron
elif [ "$3" == "3" ]
then
	 $GROMACS solvate -cp $pdb"_box".gro -cs spc216.gro -o $pdb"_solvate".gro -p topol.top
elif [ "$3" == "4" ]
then
	 $GROMACS grompp -f ions.mdp -c $pdb"_solvate".gro -p topol.top -o ions.tpr 
	#  $GROMACS dump -s ions.tpr > ions_tpr.txt
elif [ "$3" == "5" ]
then
	 printf "13" | $GROMACS genion -s ions.tpr -o $pdb"_solvate_ions".gro -p topol.top -neutral
elif [ "$3" == "6" ]
then
	 $GROMACS grompp -f em.mdp -c $pdb"_solvate_ions".gro -p topol.top -o em.tpr
	#  $GROMACS dump -s em.tpr > em_tpr.txt
elif [ "$3" == "7" ]
then
	 taskset -c 1-128 $GROMACS mdrun -v -deffnm em -nt 125
elif [ "$3" == "8" ]
then
	 $GROMACS grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr 
	#  $GROMACS dump -s nvt.tpr > nvt_tpr.txt
elif [ "$3" == "9" ]
then
	 taskset -c 1-128 $GROMACS mdrun -v -deffnm nvt -nt 125
elif [ "$3" == "10" ]
then
	 printf "16 0" | $GROMACS energy -f nvt.edr -o temperature.xvg
elif [ "$3" == "11" ]
then
	 $GROMACS grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr 
	#  $GROMACS dump -s npt.tpr > npt_tpr.txt
elif [ "$3" == "12" ]
then
	 taskset -c 1-128 $GROMACS mdrun -v -deffnm npt -nt 125
elif [ "$3" == "13" ]
then
	 printf "18 0" | $GROMACS energy -f npt.edr -o pressure.xvg
elif [ "$3" == "14" ]
then
	 printf "24 0" | $GROMACS energy -f npt.edr -o density.xvg
elif [ "$3" == "15" ]
then
	 $GROMACS grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
elif [ "$3" == "16" ]
then
	 taskset -c 1-128 $GROMACS mdrun -v -deffnm md -nt 125
elif [ "$3" == "17" ]
then
	 taskset -c 1-128 $GROMACS mdrun -v -s md.tpr -cpi md.cpt -deffnm md -append -nt 125
fi
