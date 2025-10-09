%Finally, calculate the values associated with just a stepwise model of the
%the variables in the non-model paper

Cruise=(metadata.Yearday<125)*1+(metadata.Yearday>125&metadata.Yearday<175)*2+(metadata.Yearday>175&metadata.Yearday<225)*3+(metadata.Yearday>225)*4;
Surf=metadata.Depth<2;
Bot=metadata.Depth>2;

clf
colormap default
subplot(231)
scatter(log(pprod_freference(Surf)),psbEt12(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(pprod_freference(Bot)),psbEt12(Bot),40,Cruise(Bot),'filled','d');hold on
xlabel('log primary productivity (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(a) {\itpsbE} ({\itr}="+ num2str(corr_obs(20), '%5.2f') + ", {\itp}=" + num2str(pval(20), '%5.2f') + ")", 'fontweight','bold','fontsize',14)
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(232)
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),nosZt(Surf),40,Cruise(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),nosZt(Bot),40,Cruise(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(b) {\itnosZ} ({\itr}="+ num2str(corr_obs(30), '%5.2f') + ", {\itp}=" + num2str(pval(30), '%5.2f') + ")",'fontweight','bold','fontsize',14)
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(233)
scatter(log(nitri1_freference(Surf)),HAOt(Surf),40,Cruise(Surf),'filled');hold on
scatter(log(nitri1_freference(Bot)),HAOt(Bot),40,Cruise(Bot),'filled','d');
xlabel('log ammonia oxidation (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(c) {\ithao} ({\itr}="+ num2str(corr_obs(23), '%5.2f') + ", {\itp}=" + num2str(pval(23), '%5.2f') + ")",'fontweight','bold','fontsize',14)
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(234)
scatter(log((soxo_freference(Surf)+soxno3_freference(Surf)+soxno2_freference(Surf))),dsrABoxt(Surf),40,Cruise(Surf),'filled');hold on;
scatter(log((soxo_freference(Bot)+soxno3_freference(Bot)+soxno2_freference(Bot))),dsrABoxt(Bot),40,Cruise(Bot),'filled','d');
xlabel('log hydrogen sulfide oxidation (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(d) Ox. {\itdsrAB} ({\itr}="+ num2str(corr_obs(31), '%5.2f') + ", {\itp}=" + num2str(pval(31), '%5.2f') + ")",'fontweight','bold','fontsize',14)
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(235)
scatter(log((srra_freference(Surf))'),dsrABredt(Surf),40,Cruise(Surf),'filled');hold on;
scatter(log((srra_freference(Bot))'),dsrABredt(Bot),40,Cruise(Bot),'filled','d');
xlabel('log sulfate reduction (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(e) Red. {\itdsrAB} ({\itr}="+ num2str(corr_obs(32), '%5.2f') + ", {\itp}=" + num2str(pval(32), '%5.2f') + ")",'fontweight','bold','fontsize',14)
%hc=colorbar('Location',['eastoutside']);
grid on

cmap=colormap;
colormap(cmap([1 86 171 255],:))
hc1=colorbar();
caxis([0.5 4.5]);
hc1.Position =[0.04 0.35 0.02 0.3];
xlabel(hc1,'Cruise Month');
labelsc={'Apr','Jun','Jul','Aug'}
hc1.Label.Rotation=0;
hc1.Label.Position=[0.05 0.35 0];
hc1.Direction="reverse";
hc1.Ticks=(1:4);
hc1.TickLabels=labelsc;



