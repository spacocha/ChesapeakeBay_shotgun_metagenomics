% Need to have run process_tpm and loaded haotaxat.txt
% run process_TPM_chesapeake_SPP, ches_freference_SPP, relative_rates_figure command

%Read in the taxa table
%created from the following
%in normalized_cleaned_RPKM/TPM_merged_by_KO_final_dir
%grep "Gene_locus" TPM_final_table.KO.neg.good.taxa.txt > nosZ_gammaproteobacteria.txt
%cat TPM_final_table.KO.neg.good.taxa.txt | grep "K00376" | grep "Gammaproteobacteria" >> nosZ_gammaproteobacteria.txt
%grep "Gene_locus" TPM_final_table.KO.neg.good.taxa.txt > nosZ_Bacteriodetes.txt 
%cat TPM_final_table.KO.neg.good.taxa.txt | grep "K00376" | grep "Bacteroidetes" >> nosZ_Bacteroidetes.txt
%In excel nosZ_by_taxa.xlsx, make just the headers we want and sum all
%save as just the sum for Bacteroidetes\nGammas
nosZ_taxa=readmatrix("nosZ_gamma_bacteroidetes.txt");
nosZbactert=log(nosZ_taxa(1,:))';
nosZgammat=log(nosZ_taxa(2,:)+0.01)';



%snp_hydro
clf
colormap default
subplot(221)
scatter(metadata.DO,nosZt,40,metadata.Yearday,'filled')
xlabel('log Obs. Oxygen (\muM)')
ylabel('log Abundance in nRPKM')
title('\rm(a) Obs \it{nosZ}\rm vs. Obs. DO,time')
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(222)
scatter(log(dno2_freference+dno3_freference),nosZt,40,metadata.Yearday,'filled')
xlabel('log Denitrif. rate (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(b) Obs \it{nosZ}\rm vs. Mod. Denitrif.,time')
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(223)
scatter(log(dno2_freference+dno3_freference),nosZgammat,40,metadata.Yearday,'filled')
xlabel('log Denitrif rate (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(c) Obs Gamma \it{nosZ}\rm vs. Mod. Denitrif.,time')
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(224)
scatter(log(dno2_freference+dno3_freference),nosZbactert,40,metadata.Yearday,'filled')
xlabel('log Denitrif. rate (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(d)Obs Bacteroid. \it{nosZ}\rm vs. Mod. Denitrif.,time')
%hc=colorbar('Location',['eastoutside']);
grid on

hc1=colorbar;
hc1.Position =[0.05 0.35 0.02 0.3];
ylabel(hc1,'decimal day');