#!/bin/bash -l

#SBATCH

#SBATCH --job-name=trimmo38
#SBATCH --time=2:00:00

module load trimmomatic/0.38
EX_LOC=trimmomatic
AD_LOC=/scratch/groups/sprehei1/lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
LIB_DIR=./FASTQ

ASSEM=./Trimmed_FASTQ
FWD_FQ=`awk "NR==$SLURM_ARRAY_TASK_ID" forward`
REV_FQ=`awk "NR==$SLURM_ARRAY_TASK_ID" reverse`

gunzip $LIB_DIR/${FWD_FQ}.gz
gunzip $LIB_DIR/${REV_FQ}.gz

$EX_LOC PE -threads 8 -phred33 \
$LIB_DIR/$FWD_FQ  $LIB_DIR/$REV_FQ \
$ASSEM/${FWD_FQ}_s1_pe $ASSEM/${FWD_FQ}_s1_se $ASSEM/${REV_FQ}_s2_pe $ASSEM/${REV_FQ}_s2_se \
ILLUMINACLIP:$AD_LOC:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:33 MINLEN:50

echo "Complete: trimmomatic $SLURM_ARRAY_TASK_ID"
