% Need to have run process_tpm and loaded haotaxat.txt
% run process_TPM_chesapeake_SPP, ches_freference_SPP, relative_rates_figure command

%mk_psbC
%can't subtract log values, so use these
nlpsbCt1=metabolicgenesnoatpHGnorm(konoatp==2705.1,:)';
nlpsbCt2=metabolicgenesnoatpHGnorm(konoatp==2705.2,:)';
psbCt12=log(nlpsbCt1 + nlpsbCt2);

Cruise=(metadata.Yearday<125)*1+(metadata.Yearday>125&metadata.Yearday<175)*2+(metadata.Yearday>175&metadata.Yearday<225)*3+(metadata.Yearday>225)*4;
Surf=metadata.Depth<2;
Bot=metadata.Depth>2;
%snp_hydro
clf
colormap default
subplot(231)
scatter(log(metadata.CHLA(Surf)),psbCt12(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(metadata.CHLA(Bot)),psbCt12(Bot),40,Cruise(Bot),'filled','d');
xlabel('log observed Chl a (mg m^-^3)')
ylabel('log normalized gene abundance')
title("(a) {\itpsbC} ({\itr}="+ num2str(corr_obs(36), '%5.2f') + ", {\itp}=" + num2str(pval(36), '%5.2f') + ")", 'fontweight','bold','fontsize',12)
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(232)
scatter(log(chl_freference(Surf)),psbCt12(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(chl_freference(Bot)),psbCt12(Bot),40,Cruise(Bot),'filled','d');
xlabel('log modeled Chl a (mg m^-^3)')
ylabel('log normalized gene abundance')
title("(b) {\itpsbC} ({\itr}="+ num2str(corr_obs(42), '%5.2f') + ", {\itp}=" + num2str(pval(42), '%5.2f') + ")",'fontweight','bold','fontsize',12)
%hc=colorbar('Location',['eastoutside']);
grid on
%text(-0.5,-5.9,'Circles: Surface, Diamonds: Bottom')


subplot(233)
scatter(log(pprod_freference(Surf)),psbCt12(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(pprod_freference(Bot)),psbCt12(Bot),40,Cruise(Bot),'filled','d');hold on
xlabel('log primary productivity (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(c) {\itpsbC} ({\itr}="+ num2str(corr_obs(18), '%5.2f') + ", {\itp}=" + num2str(pval(18), '%5.2f') + ")", 'fontweight','bold','fontsize',12)
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(234)
scatter(log(metadata.CHLA(Surf)),psbCt(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(metadata.CHLA(Bot)),psbCt(Bot),40,Cruise(Bot),'filled','d');
xlabel('log observed Chl a (mg m^-^3)')
ylabel('log normalized gene abundance')
title("(d) {\itCyanobacteria psbC} ({\itr}="+ num2str(corr_obs(38), '%5.2f') + ", {\itp}=" + num2str(pval(38), '%5.2f') + ")",'fontweight','bold','fontsize',12)
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(235)
scatter(log(chl_freference(Surf)),psbCt(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(chl_freference(Bot)),psbCt(Bot),40,Cruise(Bot),'filled','d');
xlabel('log modeled Chl a (mg m^-^3)')
ylabel('log normalized gene abundance')
title("(e) {\itCyanobacteria psbC} ({\itr}="+ num2str(corr_obs(39), '%5.2f') + ", {\itp}=" + num2str(pval(39), '%5.2f') + ")", 'fontweight','bold','fontsize',12)
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(236)
scatter(log(pprod_freference(Surf)),psbCt(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(pprod_freference(Bot)),psbCt(Bot),40,Cruise(Bot),'filled','d');
xlabel('log primary productivity (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(f) {\itCyanobacteria psbC} ({\itr}="+ num2str(corr_obs(40), '%5.2f') + ", {\itp}=" + num2str(pval(40), '%5.2f') + ")", 'fontweight','bold','fontsize',12)
%hc=colorbar('Location',['eastoutside']);
grid on

colormap("default")
cmap=colormap;
colormap(cmap([1 86 171 255],:))
hc1=colorbar();
caxis([0.5 4.5])
hc1.Position =[0.04 0.35 0.02 0.3];
xlabel(hc1,'Cruise Month');
labelsc={'Apr','Jun','Jul','Aug'}
hc1.Label.Rotation=0;
hc1.Label.Position=[0.05 0.35 0];
hc1.Direction="reverse";
hc1.Ticks=(1:4);
hc1.TickLabels=labelsc;