%Run process_TPM_chesapeake_SPP_rand.m
%Download these files to:
%/Users/sarahpreheim/Documents/Solexa_dir/Chesapeake_Bay/metagenomics/metagenomic_assembly/Mainstem/MATLAB
%Run the following:
%perl make_dataXY_correl_matlab2.pl lnmetabolicgenesnoatp_rand.txt konoatp_rand.txt gene_gene_correl_rand
%upload
test_uat_X=readtable("gene_metadata_corr2.X.txt");
test_uat_Y=readtable("gene_metadata_corr2.Y.txt");
dataX = table2array(test_uat_X(2:end,2:end));
dataY = table2array(test_uat_Y(2:end,2:end));
[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX',dataY', 5000, 0, .05, "rank", 1);
writematrix(pval, "pval_gene_envars2_rank.txt");
writematrix(corr_obs, "corr_gene_envars2_rank.txt");