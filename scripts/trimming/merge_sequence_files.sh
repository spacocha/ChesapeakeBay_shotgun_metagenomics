#!/bin/bash -l

#SBATCH

#SBATCH --job-name=trimmo38
#SBATCH --time=2:00:00

#merge together all of the trimmed sequence reads, whether paired or single end from both lanes of Illumina (sprehei1_146738, sprehei1_146738)
cat ../../Raw_data/sprehei1_146738/Trimmed_FASTQ/*s1_pe ../../Raw_data/sprehei1_146742/Trimmed_FASTQ/*s1_pe > all_s1_pe.fastq
cat ../../Raw_data/sprehei1_146738/Trimmed_FASTQ/*s2_pe ../../Raw_data/sprehei1_146742/Trimmed_FASTQ/*s2_pe > all_s2_pe.fastq
cat ../../Raw_data/sprehei1_146738/Trimmed_FASTQ/*se ../../Raw_data/sprehei1_146742/Trimmed_FASTQ/*se > all_se.fastq

#Use these files as input to SPADES

