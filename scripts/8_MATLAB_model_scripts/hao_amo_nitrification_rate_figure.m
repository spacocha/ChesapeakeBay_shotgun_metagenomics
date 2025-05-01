% Need to have run process_tpm and loaded haotaxat.txt
% and run command [HAOtxnum HAOtx]=sorttaxa(haotaxa,50,indenv+1)
% run process_TPM_chesapeake_SPP, ches_freference_SPP, relative_rates_figure command

%mk_amos
%can't subtract log values, so use these
nlamoAt=metabolicgenesnoatp(konoatp==10944.1,:)';
nlamoBt=metabolicgenesnoatp(konoatp==10945.1,:)';
amoABt=log(nlamoAt + nlamoBt);

%snp_hydro
clf
subplot(221)
scatter(log(metadata.NH4F),HAOt,40,metadata.Yearday,'filled')
xlabel('log Obs. Ammonium (\muM)')
ylabel('log Abundance in nRPKM')
title('\rm(a) Obs \it{hao}\rm vs. Obs Amm.,time')
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(222)
scatter(log(nitri1_freference),HAOt,40,metadata.Yearday,'filled')
xlabel('log Nitrif. rate (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(b) Obs \it{hao}\rm vs. Modeled Nitrif.,time')
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(223)
scatter(log(metadata.NH4F),amoABt,40,metadata.Yearday,'filled')
xlabel('log Obs. Ammonium (\muM)')
ylabel('log Abundance in nRPKM')
title('\rm(c) Obs \it{amoAB}\rm vs. Obs Amm,,time')
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(224)
scatter(log(nitri1_freference),amoABt,40,metadata.Yearday,'filled')
xlabel('log Nitrif. rate (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(d)Obs \it{amoAB}\rm vs. Modeled Nitrif.,time')
%hc=colorbar('Location',['eastoutside']);
grid on

hc1=colorbar;
hc1.Position =[0.05 0.35 0.02 0.3];
ylabel(hc1,'decimal day');