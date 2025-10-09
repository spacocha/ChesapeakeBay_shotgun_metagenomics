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
%divide by housekeeping genes...

Cruise=(metadata.Yearday<125)*1+(metadata.Yearday>125&metadata.Yearday<175)*2+(metadata.Yearday>175&metadata.Yearday<225)*3+(metadata.Yearday>225)*4;
Surf=metadata.Depth<2;
Bot=metadata.Depth>2;

%snp_hydro
clf
colormap default
subplot(231)
%napABt
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),napABt(Surf),40,Cruise(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),napABt(Bot),40,Cruise(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(a) {\itnapAB} ({\itr}="+ num2str(corr_obs(25), '%5.2f') + ", {\itp}=" + num2str(pval(25), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(232)
%narIt
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),narIt(Surf),40,Cruise(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),narIt(Bot),40,Cruise(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(b) {\itnarI} ({\itr}="+ num2str(corr_obs(26), '%5.2f') + ", {\itp}=" + num2str(pval(26), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(233)
%nirKt
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),nirKt(Surf),40,Cruise(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),nirKt(Bot),40,Cruise(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(c) {\itnirK} ({\itr}="+ num2str(corr_obs(27), '%5.2f') + ", {\itp}=" + num2str(pval(27), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(234)
%nirSt
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),nirSt(Surf),40,Cruise(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),nirSt(Bot),40,Cruise(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(d) {\itnirS} ({\itr}="+ num2str(corr_obs(28), '%5.2f') + ", {\itp}=" + num2str(pval(28), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

%norBCt
subplot(235)
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),norBCt(Surf),40,Cruise(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),norBCt(Bot),40,Cruise(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(e) {\itnorBC} ({\itr}="+ num2str(corr_obs(29), '%5.2f') + ", {\itp}=" + num2str(pval(29), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

%nosZt
subplot(236)
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),nosZt(Surf),40,Cruise(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),nosZt(Bot),40,Cruise(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(f) {\itnosZ} ({\itr}="+ num2str(corr_obs(30), '%5.2f') + ", {\itp}=" + num2str(pval(30), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

colormap("default")
cmap=colormap;
colormap(cmap([1 86 171 255],:))
hc1=colorbar();
caxis([0.5 4.5])
hc1.Position =[0.05 0.35 0.02 0.3];
xlabel(hc1,'Cruise Month');
labelsc={'Apr','Jun','Jul','Aug'}
hc1.Label.Rotation=0;
hc1.Label.Position=[0.08 0.35 0];
hc1.Direction="reverse";
hc1.Ticks=(1:4);
hc1.TickLabels=labelsc;