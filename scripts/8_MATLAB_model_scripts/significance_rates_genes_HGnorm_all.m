%First, make a matrix of all of the values that need to be run
%Run process_TPM_chesapeake_SPP.m first to get the genes and metadata
%Run ches_freference_SPP.m to get the modeled data
clear("pval", "corr_obs", "crit_corr", "est_alpha", "seed_state", "dataX", "dataY");


%NS
nlpetAt1=metabolicgenesnoatpHGnorm(konoatp==2634.1,:)';
nlpetAt2=metabolicgenesnoatpHGnorm(konoatp==2634.2,:)';
petAt12=log(nlpetAt1 + nlpetAt2);

dataX(1,:)=petAt12;
dataY(1,:)=log(pprod_freference);

mdlpetAt12 = fitlm(petAt12, log(pprod_freference));
results(1,1)=mdlpetAt12.Rsquared.Adjusted;

%NS
petBt=lnmetabolicgenesnoatp(konoatp==2635.1,:)';
nlpetBt1=metabolicgenesnoatpHGnorm(konoatp==2635.1,:)';
nlpetBt2=metabolicgenesnoatpHGnorm(konoatp==2635.2,:)';
petBt12=log(nlpetBt1 + nlpetBt2);

dataX(2,:)=petBt12;
dataY(2,:)=log(pprod_freference);

mdlpetBt12 = fitlm(petBt12, log(pprod_freference));
results(1,2)=mdlpetBt12.Rsquared.Adjusted;

%NS
petCt=lnmetabolicgenesnoatp(konoatp==2636.1,:)';
nlpetCt1=metabolicgenesnoatpHGnorm(konoatp==2636.1,:)';
nlpetCt2=metabolicgenesnoatpHGnorm(konoatp==2636.2,:)';
petCt12=log(nlpetCt1 + nlpetCt2);

dataX(3,:)=petCt12;
dataY(3,:)=log(pprod_freference);

mdlpetCt12 = fitlm(petCt12, log(pprod_freference));
results(1,3)=mdlpetCt12.Rsquared.Adjusted;

%pval<0.05
petDt=lnmetabolicgenesnoatp(konoatp==2637.1,:)';
nlpetDt1=metabolicgenesnoatpHGnorm(konoatp==2637.1,:)';
nlpetDt2=metabolicgenesnoatpHGnorm(konoatp==2637.2,:)';
petDt12=log(nlpetDt1 + nlpetDt2);

dataX(4,:)=petDt12;
dataY(4,:)=log(pprod_freference);

mdlpetDt12 = fitlm(petDt12, log(pprod_freference));
results(1,4)=mdlpetDt12.Rsquared.Adjusted;


%NS
petEt=lnmetabolicgenesnoatp(konoatp==2638.1,:)';
nlpetEt1=metabolicgenesnoatpHGnorm(konoatp==2638.1,:)';
nlpetEt2=metabolicgenesnoatpHGnorm(konoatp==2638.2,:)';
petEt12=log(nlpetEt1 + nlpetEt2);

dataX(5,:)=petEt12;
dataY(5,:)=log(pprod_freference);

mdlpetEt12 = fitlm(petEt12, log(pprod_freference));
results(1,5)=mdlpetEt12.Rsquared.Adjusted;

%pval<0.05
petFt=lnmetabolicgenesnoatp(konoatp==2639.1,:)';
nlpetFt1=metabolicgenesnoatpHGnorm(konoatp==2639.1,:)';
nlpetFt2=metabolicgenesnoatpHGnorm(konoatp==2639.2,:)';
petFt12=log(nlpetFt1 + nlpetFt2);

dataX(6,:)=petFt12;
dataY(6,:)=log(pprod_freference);

mdlpetFt12 = fitlm(petFt12, log(pprod_freference));
results(1,6)=mdlpetFt12.Rsquared.Adjusted;

%NS
petGt=lnmetabolicgenesnoatp(konoatp==2640.1,:)';
nlpetGt1=metabolicgenesnoatpHGnorm(konoatp==2640.1,:)';
nlpetGt2=metabolicgenesnoatpHGnorm(konoatp==2640.2,:)';
petGt12=log(nlpetGt1 + nlpetGt2);

dataX(7,:)=petGt12;
dataY(7,:)=log(pprod_freference);

mdlpetGt12 = fitlm(petGt12, log(pprod_freference));
results(1,7)=mdlpetGt12.Rsquared.Adjusted;

%NS
petHt=lnmetabolicgenesnoatp(konoatp==2641.1,:)';
nlpetHt1=metabolicgenesnoatpHGnorm(konoatp==2641.1,:)';
nlpetHt2=metabolicgenesnoatpHGnorm(konoatp==2641.2,:)';
petHt12=log(nlpetHt1 + nlpetHt2);

dataX(8,:)=petHt12;
dataY(8,:)=log(pprod_freference);

mdlpetHt12 = fitlm(petHt12, log(pprod_freference));
results(1,8)=mdlpetHt12.Rsquared.Adjusted;

petJt=lnmetabolicgenesnoatp(konoatp==8906,:)';
nlpetJt1=metabolicgenesnoatpHGnorm(konoatp==8906.1,:)';
nlpetJt2=metabolicgenesnoatpHGnorm(konoatp==8906.2,:)';
petJt12=log(nlpetJt1 + nlpetJt2);

dataX(9,:)=petJt12;
dataY(9,:)=log(pprod_freference);

mdlpetJt12 = fitlm(petJt12, log(pprod_freference));
results(1,9)=mdlpetJt12.Rsquared.Adjusted;

%NS
psaAt=lnmetabolicgenesnoatp(konoatp==2689.1,:)';
nlpsaAt1=metabolicgenesnoatpHGnorm(konoatp==2689.1,:)';
nlpsaAt2=metabolicgenesnoatpHGnorm(konoatp==2689.2,:)';
psaAt12=log(nlpsaAt1 + nlpsaAt2);

dataX(10,:)=psaAt12;
dataY(10,:)=log(pprod_freference);

mdlpsaAt12 = fitlm(psaAt12, log(pprod_freference));
results(1,10)=mdlpsaAt12.Rsquared.Adjusted;

%NS
psaBt=lnmetabolicgenesnoatp(konoatp==2690.1,:)';
nlpsaBt1=metabolicgenesnoatpHGnorm(konoatp==2690.1,:)';
nlpsaBt2=metabolicgenesnoatpHGnorm(konoatp==2690.2,:)';
psaBt12=log(nlpsaBt1 + nlpsaBt2);

dataX(11,:)=psaBt12;
dataY(11,:)=log(pprod_freference);

mdlpsaBt12 = fitlm(psaBt12, log(pprod_freference));
results(1,11)=mdlpsaBt12.Rsquared.Adjusted;

%NS
psaCt=lnmetabolicgenesnoatp(konoatp==2691.1,:)';
nlpsaCt1=metabolicgenesnoatpHGnorm(konoatp==2691.1,:)';
nlpsaCt2=metabolicgenesnoatpHGnorm(konoatp==2691.2,:)';
psaCt12=log(nlpsaCt1 + nlpsaCt2);

dataX(12,:)=psaCt12;
dataY(12,:)=log(pprod_freference);

mdlpsaCt12 = fitlm(psaCt12, log(pprod_freference));
results(1,12)=mdlpsaCt12.Rsquared.Adjusted;

%pval<0.05
psaDt=lnmetabolicgenesnoatp(konoatp==2692.1,:)';
nlpsaDt1=metabolicgenesnoatpHGnorm(konoatp==2692.1,:)';
nlpsaDt2=metabolicgenesnoatpHGnorm(konoatp==2692.2,:)';
psaDt12=log(nlpsaDt1 + nlpsaDt2);

dataX(13,:)=psaDt12;
dataY(13,:)=log(pprod_freference);

mdlpsaDt12 = fitlm(psaDt12, log(pprod_freference));
results(1,13)=mdlpsaDt12.Rsquared.Adjusted;

%pval<0.05
psaEt=lnmetabolicgenesnoatp(konoatp==2693.1,:)';
nlpsaEt1=metabolicgenesnoatpHGnorm(konoatp==2693.1,:)';
nlpsaEt2=metabolicgenesnoatpHGnorm(konoatp==2693.2,:)';
psaEt12=log(nlpsaEt1 + nlpsaEt2);

dataX(14,:)=psaEt12;
dataY(14,:)=log(pprod_freference);

mdlpsaEt12 = fitlm(psaEt12, log(pprod_freference));
results(1,14)=mdlpsaEt12.Rsquared.Adjusted;

%NS
psaFt=lnmetabolicgenesnoatp(konoatp==2694.1,:)';
nlpsaFt1=metabolicgenesnoatpHGnorm(konoatp==2694.1,:)';
nlpsaFt2=metabolicgenesnoatpHGnorm(konoatp==2694.2,:)';
psaFt12=log(nlpsaFt1 + nlpsaFt2);

dataX(15,:)=psaFt12;
dataY(15,:)=log(pprod_freference);

mdlpsaFt12 = fitlm(psaFt12, log(pprod_freference));
results(1,15)=mdlpsaFt12.Rsquared.Adjusted;

%NS
psbAt=lnmetabolicgenesnoatp(konoatp==2703.1,:)';
nlpsbAt1=metabolicgenesnoatpHGnorm(konoatp==2703.1,:)';
nlpsbAt2=metabolicgenesnoatpHGnorm(konoatp==2703.2,:)';
psbAt12=log(nlpsbAt1 + nlpsbAt2);

dataX(16,:)=psbAt12;
dataY(16,:)=log(pprod_freference);

mdlpsbAt12 = fitlm(psbAt12, log(pprod_freference));
results(1,16)=mdlpsbAt12.Rsquared.Adjusted;

%pval<0.05
psbBt=lnmetabolicgenesnoatp(konoatp==2704.1,:)';
nlpsbBt1=metabolicgenesnoatpHGnorm(konoatp==2704.1,:)';
nlpsbBt2=metabolicgenesnoatpHGnorm(konoatp==2704.2,:)';
psbBt12=log(nlpsbBt1 + nlpsbBt2);

dataX(17,:)=psbBt12;
dataY(17,:)=log(pprod_freference);

mdlpsbBt12 = fitlm(psbBt12, log(pprod_freference));
results(1,17)=mdlpsbBt12.Rsquared.Adjusted;

%pval<0.05
psbCt=lnmetabolicgenesnoatp(konoatp==2705.1,:)';
nlpsbCt1=metabolicgenesnoatpHGnorm(konoatp==2705.1,:)';
nlpsbCt2=metabolicgenesnoatpHGnorm(konoatp==2705.2,:)';
psbCt12=log(nlpsbCt1 + nlpsbCt2);

dataX(18,:)=psbCt12;
dataY(18,:)=log(pprod_freference);

mdlpsbCt12 = fitlm(psbCt12, log(pprod_freference));
results(1,18)=mdlpsbCt12.Rsquared.Adjusted;

%NS
psbDt=lnmetabolicgenesnoatp(konoatp==2706.1,:)';
nlpsbDt1=metabolicgenesnoatpHGnorm(konoatp==2706.1,:)';
nlpsbDt2=metabolicgenesnoatpHGnorm(konoatp==2706.2,:)';
psbDt12=log(nlpsbDt1 + nlpsbDt2);

dataX(19,:)=psbDt12;
dataY(19,:)=log(pprod_freference);

mdlpsbDt12 = fitlm(psbDt12, log(pprod_freference));
results(1,19)=mdlpsbDt12.Rsquared.Adjusted;

%pval<0.05
psbEt=lnmetabolicgenesnoatp(konoatp==2707.1,:)';
nlpsbEt1=metabolicgenesnoatpHGnorm(konoatp==2707.1,:)';
nlpsbEt2=metabolicgenesnoatpHGnorm(konoatp==2707.2,:)';
psbEt12=log(nlpsbEt1 + nlpsbEt2);

dataX(20,:)=psbEt12;
dataY(20,:)=log(pprod_freference);

mdlpsbEt12 = fitlm(psbEt12, log(pprod_freference));
results(1,20)=mdlpsbEt12.Rsquared.Adjusted;

%NS
psbFt=lnmetabolicgenesnoatp(konoatp==2708.1,:)';
nlpsbFt1=metabolicgenesnoatpHGnorm(konoatp==2708.1,:)';
nlpsbFt2=metabolicgenesnoatpHGnorm(konoatp==2708.2,:)';
psbFt12=log(nlpsbFt1 + nlpsbFt2);

dataX(21,:)=psbFt12;
dataY(21,:)=log(pprod_freference);

mdlpsbFt12 = fitlm(psbFt12, log(pprod_freference));
results(1,21)=mdlpsbFt12.Rsquared.Adjusted;

%pval<0.05
psbOt=lnmetabolicgenesnoatp(konoatp==2716.1,:)';
nlpsbOt1=metabolicgenesnoatpHGnorm(konoatp==2716.1,:)';
nlpsbOt2=metabolicgenesnoatpHGnorm(konoatp==2716.2,:)';
psbOt12=log(nlpsbOt1 + nlpsbOt2);

dataX(22,:)=psbOt12;
dataY(22,:)=log(pprod_freference);

mdlpsbOt12 = fitlm(psbOt12, log(pprod_freference));
results(1,22)=mdlpsbOt12.Rsquared.Adjusted;

%pval<0.05
dataX(23,:)=HAOt;
dataY(23,:)=log(nitri1_freference);
%pval 0
%corr 0.7137

mdlHAOt = fitlm(HAOt, log(nitri1_freference));
results(1,23)=mdlHAOt.Rsquared.Adjusted;

%NS
nlamoAt=metabolicgenesnoatpHGnorm(konoatp==10944.1,:)';
nlamoBt=metabolicgenesnoatpHGnorm(konoatp==10945.1,:)';
amoABt=log(nlamoAt + nlamoBt);

dataX(24,:)=amoABt;
dataY(24,:)=log(nitri1_freference);

mdlamoABt = fitlm(amoABt, log(nitri1_freference));
results(1,24)=mdlamoABt.Rsquared.Adjusted;
%pval 0.5672
%corr 0.2909

%pval<0.05
nlnapAt=metabolicgenesnoatpHGnorm(konoatp==2567,:)';
nlnapBt=metabolicgenesnoatpHGnorm(konoatp==2568,:)';
napABt=log(nlnapAt + nlnapBt);

dataX(25,:)=napABt;
dataY(25,:)=log((dno3_freference+dno2_freference));

mdlnapABt = fitlm(napABt, log((dno3_freference+dno2_freference)));
results(1,25)=mdlnapABt.Rsquared.Adjusted;


%pval<0.05
narIt=lnmetabolicgenesnoatp(konoatp==374,:)';
dataX(26,:)=narIt;
dataY(26,:)=log((dno3_freference+dno2_freference));

mdlnarIt = fitlm(narIt, log((dno3_freference+dno2_freference)));
results(1,26)=mdlnarIt.Rsquared.Adjusted;

%NS
dataX(27,:)=nirKt;
dataY(27,:)=log((dno3_freference+dno2_freference));

mdlnirKt = fitlm(nirKt, log((dno3_freference+dno2_freference)));
results(1,27)=mdlnirKt.Rsquared.Adjusted;

%NS
dataX(28,:)=nirSt;
dataY(28,:)=log((dno3_freference+dno2_freference));

mdlnirSt = fitlm(nirSt, log((dno3_freference+dno2_freference)));
results(1,28)=mdlnirSt.Rsquared.Adjusted;

%pval<0.05
nlnorBt=metabolicgenesnoatpHGnorm(konoatp==4561,:)';
nlnorCt=metabolicgenesnoatpHGnorm(konoatp==2305,:)';
norBCt=log(nlnorBt + nlnorCt);

dataX(29,:)=norBCt;
dataY(29,:)=log((dno3_freference+dno2_freference));

mdlnorBCt = fitlm(norBCt, log((dno3_freference+dno2_freference)));
results(1,29)=mdlnorBCt.Rsquared.Adjusted;

%pval<0.05
dataX(30,:)=nosZt;
dataY(30,:)=log((dno3_freference+dno2_freference));

mdlnosZt = fitlm(nosZt, log((dno3_freference+dno2_freference)));
results(1,30)=mdlnosZt.Rsquared.Adjusted;
%pval 0
%corr 0.7208


%pval<0.05
nldsrAoxt=metabolicgenesnoatpHGnorm(konoatp==11180.1,:)';
nldsrBoxt=metabolicgenesnoatpHGnorm(konoatp==11181.1,:)';
dsrABoxt=log(nldsrAoxt+ nldsrBoxt);
dataX(31,:)=dsrABoxt;
dataY(31,:)=log((soxo_freference + soxno2_freference + soxno3_freference));

mdldsrABoxt = fitlm(dsrABoxt, log((soxo_freference + soxno2_freference + soxno3_freference)));
results(1,31)=mdldsrABoxt.Rsquared.Adjusted;
%pval 0.0072
%corr 0.5450

%NS
%nlsoxAt=metabolicgenesnoatpHGnorm(konoatp==17222,:)';
%nlsoxXt=metabolicgenesnoatpHGnorm(konoatp==17223,:)';
%soxAXt=log(nlsoxAt+ nlsoxXt);
%dataX(32,:)=soxAXt;
%dataY(32,:)=log((soxo_freference + soxno2_freference + soxno3_freference));

%mdlsoxAXt = fitlm(soxAXt, log((soxo_freference + soxno2_freference + soxno3_freference)));
%results(1,32)=mdlsoxAXt.Rsquared.Adjusted;

%pval<0.05 (opposite direction)
%nlsoxCt=metabolicgenesnoatpHGnorm(konoatp==17225,:)';
%nlsoxDt=metabolicgenesnoatpHGnorm(konoatp==22622,:)';
%soxCDt=log(nlsoxCt+ nlsoxDt);
%dataX(33,:)=soxCDt;
%dataY(33,:)=log((soxo_freference + soxno2_freference + soxno3_freference));

%mdlsoxCDt = fitlm(soxCDt, log((soxo_freference + soxno2_freference + soxno3_freference)));
%results(1,33)=mdlsoxCDt.Rsquared.Adjusted;

%NS
%nlsoxBt=metabolicgenesnoatpHGnorm(konoatp==17224,:)';
%dataX(34,:)=soxBt;
%dataY(34,:)=log((soxo_freference + soxno2_freference + soxno3_freference));

%mdlsoxBt = fitlm(soxBt, log((soxo_freference + soxno2_freference + soxno3_freference)));
%results(1,34)=mdlsoxBt.Rsquared.Adjusted;

%NS
%nlsoxYt=metabolicgenesnoatpHGnorm(konoatp==17226,:)';
%dataX(35,:)=soxYt;
%dataY(35,:)=log((soxo_freference + soxno2_freference + soxno3_freference));

%mdlsoxYt = fitlm(soxYt, log((soxo_freference + soxno2_freference + soxno3_freference)));
%results(1,35)=mdlsoxYt.Rsquared.Adjusted;

%NS
%nlsoxZt=metabolicgenesnoatpHGnorm(konoatp==17227,:)';
%dataX(36,:)=soxZt;
%dataY(36,:)=log((soxo_freference + soxno2_freference + soxno3_freference));

%mdlsoxZt = fitlm(soxZt, log((soxo_freference + soxno2_freference + soxno3_freference)));
%results(1,36)=mdlsoxZt.Rsquared.Adjusted;

%NS
%mk_dsrred
%can't subtract log values, so use these
nldsrAredt=metabolicgenesnoatpHGnorm(konoatp==11180.2,:)';
nldsrBredt=metabolicgenesnoatpHGnorm(konoatp==11181.2,:)';
dsrABredt=log(nldsrAredt + nldsrBredt);

%NS
dataX(32,:)=dsrABredt;
dataY(32,:)=log(srra_freference');

mdldsrABredt = fitlm(dsrABredt, log(srra_freference'));
results(1,32)=mdldsrABredt.Rsquared.Adjusted;
%pval 0.6432
%corr 0.2757

dataX(33,:)=nosZt;
dataY(33,:)=log(metadata.DO);

dataX(34,:)=nosZbactert;
dataY(34,:)=log((dno3_freference+dno2_freference));


dataX(35,:)=nosZgammat;
dataY(35,:)=log((dno3_freference+dno2_freference));

dataX(36,:)=psbCt12;
dataY(36,:)=log(metadata.CHLA);

dataX(37,:)=psbCt12;
dataY(37,:)=log(chl_freference);

%cyanobacteria
dataX(38,:)=psbCt;
dataY(38,:)=log(metadata.CHLA);

dataX(39,:)=psbCt;
dataY(39,:)=log(chl_freference);

dataX(40,:)=psbCt;
dataY(40,:)=log(pprod_freference);

dataX(41,:)=HAOt;
dataY(41,:)=log(metadata.NH4F);


dataX(42,:)=amoABt;
dataY(42,:)=log(metadata.NH4F);


dataX(43,:)=log(metadata.DO);
dataY(43,:)=log(dno3_freference+dno2_freference);


[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX',dataY');
writematrix(pval, "pval_rates_genes_HGnorm108_all_plus.txt");
writematrix(corr_obs, "corr_rates_genes_HGnorm108_all_plus.txt");
writematrix(results, "results_R2_rates_genes_HGnorm108_all_plus.txt");

