%Finally, calculate the values associated with just a stepwise model of the
%the variables in the non-model paper

Cruise=(metadata.Yearday<125)*1+(metadata.Yearday>125&metadata.Yearday<175)*2+(metadata.Yearday>175&metadata.Yearday<225)*3+(metadata.Yearday>225)*4;
Surf=metadata.Depth<2;
Bot=metadata.Depth>2;

clf
colormap("copper")
subplot(231)
scatter(log(pprod_freference(Surf)),psbEt12(Surf),40,metadata.Lat(Surf),'filled');hold on
scatter(log(pprod_freference(Bot)),psbEt12(Bot),40,metadata.Lat(Bot),'filled','d');hold on
xlabel('log primary productivity (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(a) {\itpsbE} (\itr="+ num2str(corr_obs(20), '%5.3f') + " {\itp}=" + num2str(pval(20), '%5.3f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(232)
scatter(log((dno2_freference(Surf)+dno3_freference(Surf))),nosZt(Surf),40,metadata.Lat(Surf),'filled'); hold on;
scatter(log((dno2_freference(Bot)+dno3_freference(Bot))),nosZt(Bot),40,metadata.Lat(Bot),'filled','d')
xlabel('log denitrification (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(b) {\itnosZ} ({\itr}="+ num2str(corr_obs(30), '%5.3f') + ", {\itp}=" + num2str(pval(30), '%5.3f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(233)
scatter(log(nitri1_freference(Surf)),HAOt(Surf),40,metadata.Lat(Surf),'filled');hold on
scatter(log(nitri1_freference(Bot)),HAOt(Bot),40,metadata.Lat(Bot),'filled','d');
xlabel('log nitrification (mmol m^{-3} day^{-1})')
ylabel('log relative gene abundance')
title("(c) {\ithao} (\itr="+ num2str(corr_obs(23), '%5.3f') + " {\itp}=" + num2str(pval(23), '%5.3f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(234)
scatter(log((soxo_freference(Surf)+soxno3_freference(Surf)+soxno2_freference(Surf))),dsrABoxt(Surf),40,metadata.Lat(Surf),'filled');hold on;
scatter(log((soxo_freference(Bot)+soxno3_freference(Bot)+soxno2_freference(Bot))),dsrABoxt(Bot),40,metadata.Lat(Bot),'filled','d');
xlabel('log hydrogen sulfide oxidation (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(d) Ox. {\itdsrAB} (\itr="+ num2str(corr_obs(31), '%5.3f') + " {\itp}=" + num2str(pval(31), '%5.3f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

subplot(235)
scatter(log((srra_freference(Surf))'),dsrABredt(Surf),40,metadata.Lat(Surf),'filled');hold on;
scatter(log((srra_freference(Bot))'),dsrABredt(Bot),40,metadata.Lat(Bot),'filled','d');
xlabel('log sulfate reduction (mmol m^{-3} day^{-1})')
ylabel('log normalized gene abundance')
title("(e) Red. {\itdsrAB} (\itr="+ num2str(corr_obs(32), '%5.3f') + " {\itp}=" + num2str(pval(32), '%5.3f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on

hc1=colorbar;
hc1.Position =[0.05 0.35 0.02 0.3];
ylabel(hc1,'Latitude');



