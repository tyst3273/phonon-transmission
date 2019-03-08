#!/bin/bash
#PBS -N lammps
#PBS -l nodes=1:ppn=32
#PBS -l walltime=24:00:00
###PBS -l mem=256m
#PBS -e lammps.err
#PBS -o lammps.log

#mypath=$(pwd | head)
#cd $mypath
cd "$PBS_O_WORKDIR"
NPROCS=`wc -l < $PBS_NODEFILE`

mpirun -np $NPROCS /mnt/raid5/ty/program_files/lammps-12Dec18/src/lmp_ty -var SEED $RANDOM -v langSEED $RANDOM -in in.* 
