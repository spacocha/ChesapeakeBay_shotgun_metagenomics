#Divide dsrAB into oxidative and reductive types based on previous analysis
#Download Supplemental File S2 from
#Database of oxidative and reductive types was taken from Supplemental File S2 of
#Phylogenetic and environmental diversity of DsrAB-type dissimilatory (bi) sulfite reductases
#
#A. L. Muller, K. U. Kjeldsen, T. Rattei, M. Pester and A. Loy
#
#ISME Journal 2015 Vol. 9 Issue 5 Pages 1152-1165
#
#Accession Number: WOS:000353354100009 DOI: 10.1038/ismej.2014.208
#
#https://static-content-springer-com.proxy1.library.jhu.edu/esm/art%3A10.1038%2Fismej.2014.208/MediaObjects/41396_2015_BFismej2014208_MOESM45_ESM.txt
#This is an aligned file dsrab_dome_nuc.fasta
#Use BLAST to identify which type
ml blast/2.13.0

#Get all of the K11180 genes
cat  ../GhostKhoala_redo/*[0-9].ko.txt | grep  "K11180" > K11180_list.txt
#Make those into a fasta file
perl tab2fasta_from_list.pl ../assembledContigs.all_genes.redo.fa K11180_list.txt 1 > assembledContigs.all_genes.redo.K11180.fasta

#Remove the alignment and names
perl fix_fasta_redox_type.pl dsrab_dome_nuc_remove.fasta > dsrab_dome_nuc_remove_rename.fasta

#Prepare Blast database from the SI S2 fixed file
makeblastdb -dbtype nucl -in dsrab_dome_nuc_remove_rename.fasta

blastn -db dsrab_dome_nuc_remove_rename.fasta -query assembledContigs.all_genes.redo.K11180.fasta -outfmt 6 -max_target_seqs 1 -out assembledContigs.all_genes.redo.K11180.blast

#get the relative abundance and add info to the table
perl divide_merged_KO_numbers_blast.pl ../normalized_cleaned_RPKM/TPM_merged_by_KO_final_dir/TPM_final_table.KO.neg.good.taxa.txt K11180 assembledContigs.all_genes.redo.K11180.blast.ox assembledContigs.all_genes.redo.K11180.blast.red > TPM_final_table.KO.neg.merged.divided.K11180.txt

perl divide_merged_KO_numbers_blast.pl ../normalized_cleaned_RPKM/TPM_merged_by_KO_final_dir/TPM_final_table.KO.neg.good.taxa.txt K11181 assembledContigs.all_genes.redo.K11181.blast.ox assembledContigs.all_genes.redo.K11181.blast.red > TPM_final_table.KO.neg.merged.divided.K11181.txt

#key
#K11180.1 (dsrA oxidative)
#K11180.2 (dsrA reductive)
#K11181.1 (dsrB oxidative)
#K11181.2 (dsrB reductive)


