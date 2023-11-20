#Working in MATLAB
#23.2.0.2436196 (R2023b) Update 4
#web interface at https://matlab.mathworks.com/
#JHU License

#Process the metagenomic data
#from TPM_final_table.KO.neg.merged.divided.final.txt 
process_TPM_chesapeake_SPP.m

#Output is
#writematrix(lnmetabolicgenesnoatp, "lnmetabolicgenesnoatp.txt"); (log OTU table, without housekeeping genes)
#writematrix(konoatp',"konoatp.txt"); (KO numbers for each, without housekeeping)
#writematrix(uat, "uat.txt"); (metabolic gene EOFs)
#writematrix(uct, "uct.txt"); (Common genes EOF)
#writematrix(lnhousekeepinggenes,"lnhousekeepinggenes.txt"); (log OTU table with only housekeeping genes)
#writematrix(kohousekeeping,"kohousekeeping.txt"); (Non-log housekeeping genes)
#writematrix(lnmetabolicgenesnoatpanom, "lnmetabolicgenesnoatpanom.txt"); #normalized log-transformed OTU table)

#Run the regression modeling
#Plots are generated, but not used in the paper
regression_models_merged.m

#Make the EOF figures used in Part I
mk_EOFs.m

#run the statistical tests
#with multiple hypothesis testing corrections
#All based on the previously developed code
mult_comp_perm_corr.m

#Gene vs. gene correlations
#to properly format dataX and dataY tables
perl make_dataXY_correl_matlab3.pl lnmetabolicgenesnoatp.txt konoatp.txt exclude_list_gene_gene_correlation.txt2 gene_gene_correl

#Run subunit and gene-gene  correlation analysis with mult_comp_perm_corr.m for all pairs of genes
#Some genes are excluded (based on list)
gene_genes_correlations.m

#Parse results adding the gene names to the output
perl parse_results2.pl corr_gene_gene.txt gene_gene_correl.log.txt > corr_gene_gene.results.txt
perl parse_results2.pl pval_gene_gene.txt gene_gene_correl.log.txt > pval_gene_gene.results.txt

#Identify which environmental variables are correlated with EOFs
perl make_dataXY_correl_matlab_uat_envars.pl uat.txt MATLAB_metadata_3.txt EOF_metadata_correl

#Run mult_comp_perm_corr.m with this output
EOF_envars_correlations.m

#Parse results adding the gene names to the output
perl parse_results2.pl corr_EOF_metadata.txt EOF_metadata_correl.log.txt > corr_EOF_metadata.results.txt
perl parse_results2.pl pval_EOF_metadata.txt EOF_metadata_correl.log.txt > pval_EOF_metadata.results.txt

#Identify which genes are correlated with EOFs
#Set up the gene data from MATLAB output
perl make_dataXY_correl_matlab_uat_genes.pl uat.txt lnmetabolicgenesnoatp.txt konoatp.txt test_uat

#Run mult_comp_perm_corr.m with this output
EOF_genes_correlations.m

#parse results
perl parse_results2.pl corr_EOF_gene.txt EOF_gene_correl.log.txt > corr_EOF_gene.results.txt
perl parse_results2.pl pval_EOF_gene.txt EOF_gene_correl.log.txt > pval_EOF_gene.results.txt

#Also modeling
process_TPMneg_chesapeake.m

#This requires
#pco2_calc_b.m
#This file also runs and makes the model envar observations
#ches_freference.m

#This is some data that also needs to be loaded
#TPM_final_negrem_1.mat

#Analyze the transcriptomics data
#Makes the transcript figure
process_hewson.m

#This is also the transcript data that needs to be loaded
#HewsonTPMfinaltable.mat


#Output is the Chl a figure, but other figures would be similar


