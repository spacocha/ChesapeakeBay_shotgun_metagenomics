%Run process_TPM_chesapeake_SPP_rand.m
%Download these files to:
%/Users/sarahpreheim/Documents/Solexa_dir/Chesapeake_Bay/metagenomics/metagenomic_assembly/Mainstem/MATLAB
%Run the following:
%perl make_dataXY_correl_matlab_uat_genes.pl uat.txt lnmetabolicgenesnoatp.txt konoatp.txt uat_gene
%upload
test_uat_X=readtable("uat_gene.X.txt");
test_uat_Y=readtable("uat_gene.Y.txt");
dataX = table2array(test_uat_X(2:end,2:end));
dataY = table2array(test_uat_Y(2:end,2:end));
[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX',dataY',5000, 0, .05, "rank", 1);
writematrix(pval, "pval_uat_gene_rank.txt");
writematrix(corr_obs, "corr_uat_gene_rank.txt");