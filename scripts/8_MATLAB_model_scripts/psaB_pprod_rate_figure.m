% Need to have run process_tpm and loaded haotaxat.txt
% run process_TPM_chesapeake_SPP, ches_freference_SPP, relative_rates_figure command

%mk_psaB
%can't subtract log values, so use these
nlpsaBt1=metabolicgenesnoatp(konoatp==2690.1,:)';
nlpsaBt2=metabolicgenesnoatp(konoatp==2690.2,:)';
psaBt12=log(nlpsaBt1 + nlpsaBt2);

%snp_hydro
clf
colormap default
subplot(231)
scatter(log(metadata.CHLA),psaBt12,40,metadata.Yearday,'filled')
xlabel('log Chl (mg m^-^3)')
ylabel('log Abundance in nRPKM')
title('\rm(a)Obs \it{psaB}\rm vs. Obs Chl, time')
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(232)
scatter(log(chl_freference),psaBt12,40,metadata.Yearday,'filled')
xlabel('log Chl (mg m^-^3)')
ylabel('log Abundance in nRPKM')
title('\rm(b)Obs \it{psaB}\rm vs. Mod. Chl,time')
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(233)
scatter(log(pprod_freference),psaBt12,40,metadata.Yearday,'filled')
xlabel('log Prim. Prod. (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(c)Obs \it{psaB}\rm vs. Mod. Prim. Prod.,time')
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(234)
scatter(log(metadata.CHLA),psaBt,40,metadata.Yearday,'filled')
xlabel('log Chl (mg m^-^3)')
ylabel('log Abundance in nRPKM')
title('\rm(d)Obs \it{psaB}\rm vs. Obs. Chl.,time')
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(235)
scatter(log(chl_freference),psaBt,40,metadata.Yearday,'filled')
xlabel('log Chl (mg m^-^3)')
ylabel('log Abundance in nRPKM')
title('\rm(e)Obs \it{psaB}\rm vs. Mod. Chl.,time')
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(236)
scatter(log(pprod_freference),psaBt,40,metadata.Yearday,'filled')
xlabel('log Prim. Prod (mmol m^-^3 d^-^1)')
ylabel('log Abundance in nRPKM')
title('\rm(f)Obs \it{psaB}\rm vs. Mod. Prim. Prod,time')
%hc=colorbar('Location',['eastoutside']);
grid on

hc1=colorbar;
hc1.Position =[0.05 0.35 0.02 0.3];
ylabel(hc1,'decimal day');