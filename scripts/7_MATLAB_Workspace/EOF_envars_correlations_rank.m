%Run process_TPM_chesapeake_SPP.m
%Download these files to:
%/Users/sarahpreheim/Documents/Solexa_dir/Chesapeake_Bay/metagenomics/metagenomic_assembly/Mainstem/MATLAB
%Run the following:
%perl make_dataXY_correl_matlab3.pl uat.txt MATLAB_metadata_3.txt EOF_envar_correl
%upload
test_uat_X=readtable("EOF_metadata_correl.X.txt");
test_uat_Y=readtable("EOF_metadata_correl.Y.txt");
dataX = table2array(test_uat_X(2:end,2:end));
dataY = table2array(test_uat_Y(2:end,2:end));
[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX',dataY', 5000, 0, .05, "rank", 1);
writematrix(pval, "pval_EOF_metadata_rank.txt");
writematrix(corr_obs, "corr_EOF_metadata_rank.txt");