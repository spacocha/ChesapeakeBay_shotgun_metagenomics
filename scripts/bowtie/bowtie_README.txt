This document describes how to run bowtie to map reads to contigs.
Notes
This method of analysis is best when you have run trimmomatic to remove adapters and to trim your reads before assembly. Bowtie can take both paired and single end data simultaneously, which is something that other programs, like salmon, doesn't do well.

This does not map reads as well as salmon (salmon is more sensitive and seems to map more reads), but in an analysis of how data changes over our entire dataset from both salmon and bowtie mapping, the results were similar. However, when comparing the results for a small subset of reads, salmon mapped more reads (unknown whether either was mapped correctly).

Overview
Take reads in fastq format and assembly in fasta format, map reads to assembly with bowtie, then use samtools and bedtools to identify how the reads map to gene locations, predicted by prodigal. These scrips also assume you have KEGG IDs for the genes, which can also be mapped.

Conditions
1.) Assembled contigs in fasta format
2.) Both single end and paired end reads (e.g. from trimmomatic)
3.) Gene calls from Prodigal
4.) KEGG IDs and taxonomic assignments for genes from BlastKoala

Steps
1.) Use bowtie.sh script to make the build (i.e. bowtie2-build) from assembled contigs
2.) Prepare lists of the forward reads (in this case which have format *1.fastq forward, *2..fastq reverse) as follows
ls <path to forward reads>/*1.fastq > forward
ls <path to reverse reads>/*2.fastq > reverse
ls <path to reverse reads>/*s1_se.fastq > se1
ls <path to reverse reads>/*s2_se.fastq > se2
However it is done, you need to make sure they have the same order for each sequenced file (which might be an issue for a small sample or subset, which might be missing one of these files)

3.) Make a list of base names somehow. In this case, taking one of the files and removing the beginning and end parts
perl ~/bin/beagle_bin/find_remove.pl forward ../../../salmon_key_KO/Trimmed_FASTQ/ > output
perl ~/bin/beagle_bin/find_remove.pl output _1.fastq_s1_pe_1.fastq > output2
mv output2 output

4.) Run the bowtie mapping program as an array on MARCC (this will have to be configured differently for different systems or if not using SLURM or arrays)
sbatch  -o slurm-%A_%a.out --array=1-47 bowtie2_array_5.sh

5.) Merge the output, making sure to use extra memory
interact -t 2:00:00 -p shared -c 24
ls *q20.TPM.txt > TPM_list.txt 
merge_quant_TPM.pl TPM_list.txt > TPM_final_table.txt
ls ../../../Prodigal_gene_calls/GhostKhoala_redo/*.ko.txt > KOlist
perl add_KEGG_to_table.pl KOlist TPM_final_table.txt > TPM_final_table.KO.txt
ls ../../../Prodigal_gene_calls/GhostKhoala_redo/*.taxa* > taxalist
perl add_taxa_to_table.pl taxalist TPM_final_table.txt > TPM_final_table.KO.taxa.txt
perl merge_by_KO.pl TPM_final_table.KO.txt > TPM_merged_by_KO_final.txt


