#!/bin/bash -l

#SBATCH

#SBATCH --job-name=SPDMainstem
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --partition=lrgmem

SPADES=~/work/lib/SPAdes-3.13.1-Linux/bin/spades.py
FWD_FQ=./all_s1_pe.fastq
REV_FQ=./all_s2_pe.fastq
SE_FQ=./all_se.fastq
OUT_DIR=assembledContigs_redo

free -mh
getconf _NPROCESSORS_ONLN

$SPADES --meta -1 $FWD_FQ -2 $REV_FQ -s $SE_FQ -o $OUT_DIR --threads 48 --memory 950

