% Need to have run process_tpm and loaded haotaxat.txt
% run process_TPM_chesapeake_SPP, ches_freference_SPP, relative_rates_figure command

%mk_adsrox
%can't subtract log values, so use these
nldsrAoxt=metabolicgenesnoatp(konoatp==11180.1,:)';
nldsrBoxt=metabolicgenesnoatp(konoatp==11181.1,:)';
dsrABoxt=log(nldsrAoxt + nldsrBoxt);

%mk_dsrred
%can't subtract log values, so use these
nldsrAredt=metabolicgenesnoatp(konoatp==11180.2,:)';
nldsrBredt=metabolicgenesnoatp(konoatp==11181.2,:)';
dsrABredt=log(nldsrAredt + nldsrBredt);

%snp_hydro
clf
colormap hot
subplot(221)
scatter(log(srra_freference),log(h2s_freference),40,DO_freference,'filled')
xlabel('log Mod. Sulf. Red. (mmol m^-^3 d^-^1)')
ylabel('log Mod. Sulfide (\muM)')
title('\rm(a) Mod. Sulfide vs. Sulf. Red.,DO')
colororder("reef")
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(222)
scatter(log(soxo_freference+soxno3_freference+soxno2_freference),log(h2s_freference),40,DO_freference,'filled')
xlabel('log Mod. Sulf. Ox. (mmol m^-^3 d^-^1)')
ylabel('log Mod. Sulfide (\muM)')
title('\rm(b) Mod. Sulfide vs. Sulf. Ox.,DO')
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(223)
scatter(log(srra_freference),dsrABredt,40,DO_freference,'filled')
xlabel('log Mod. Sulf. Red. (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(c) Obs Red. \it{dsrAB}\rm vs. Mod. Sulf. Red.,DO')
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(224)
scatter(log(soxo_freference+soxno3_freference+soxno2_freference),dsrABoxt,40,DO_freference,'filled')
xlabel('log Mod. Sulf Ox. (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(d) Obs Ox. \it{dsrAB}\rm vs. Mod. Sulf. Ox.,DO')
%hc=colorbar('Location',['eastoutside']);
grid on

hc1=colorbar;
hc1.Position =[0.05 0.35 0.02 0.3];
ylabel(hc1,'oxygen (\mum)');