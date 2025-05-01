%First, make a matrix of all of the values that need to be run
%Run process_TPM_chesapeake_SPP.m first to get the genes and metadata
%Run ches_freference_SPP.m to get the modeled data
clear("pval", "corr_obs", "crit_corr", "est_alpha", "seed_state");

%mk_psaB
%can't subtract log values, so use these
nlpsaBt1=metabolicgenesnoatp(konoatp==2690.1,:)';
nlpsaBt2=metabolicgenesnoatp(konoatp==2690.2,:)';
psaBt12=log(nlpsaBt1 + nlpsaBt2);

%Begin to assemble the two matrices with X and Y values to compare
%in separate matrices
%Obs, modeled, rate
dataX(1,:)=psaBt12;
dataY(1,:)=log(metadata.CHLA);
%pval 0
%corr 0.8551

dataX(2,:)=psaBt12;
dataY(2,:)=log(chl_freference);
%pval 0.1624
%corr 0.3859


dataX(3,:)=psaBt12;
dataY(3,:)=log(pprod_freference);
%pval 0.0076
%corr 0.5347

dataX(4,:)=psaBt;
dataY(4,:)=log(metadata.CHLA);
%pval 0.5588
%corr 0.2927

dataX(5,:)=psaBt;
dataY(5,:)=log(chl_freference);
%pval 0
%corr 0.6542


dataX(6,:)=psaBt;
dataY(6,:)=log(pprod_freference);
%pval 0.0064
%corr 0.5498

dataX(7,:)=HAOt;
dataY(7,:)=log(metadata.NH4F);
%pval 0
%corr 0.7748

dataX(8,:)=HAOt;
dataY(8,:)=log(nh4_freference);
%pval 0
%corr 0.6797

dataX(9,:)=HAOt;
dataY(9,:)=log(nitri1_freference);
%pval 0
%corr 0.7137

nlamoAt=metabolicgenesnoatp(konoatp==10944.1,:)';
nlamoBt=metabolicgenesnoatp(konoatp==10945.1,:)';
amoABt=log(nlamoAt + nlamoBt);

dataX(10,:)=amoABt;
dataY(10,:)=log(metadata.NH4F);
%pval 0.9160
%corr 0.2091


dataX(11,:)=amoABt;
dataY(11,:)=log(nh4_freference);
%pval 0.6500
%corr  0.2742

dataX(12,:)=amoABt;
dataY(12,:)=log(nitri1_freference);
%pval 0.5672
%corr 0.2909

dataX(13,:)=nosZt;
dataY(13,:)=log(metadata.DO);
%pval 0
%corr -0.6177

dataX(14,:)=nosZt;
dataY(14,:)=log(DO_freference);
%pval 8.0000e-04
%corr -0.5936


dataX(15,:)=nosZt;
dataY(15,:)=log((dno3_freference+dno2_freference));
%pval 0
%corr 0.7208

%mk_adsrox
%can't subtract log values, so use these
nldsrAoxt=metabolicgenesnoatp(konoatp==11180.1,:)';
nldsrBoxt=metabolicgenesnoatp(konoatp==11181.1,:)';
dsrABoxt=log(nldsrAoxt + nldsrBoxt);

dataX(16,:)=dsrABoxt;
dataY(16,:)=log(h2s_freference);
%pval 0.0068
%corr 0.5487

dataX(17,:)=dsrABoxt;
dataY(17,:)=log((soxo_freference + soxno2_freference + soxno3_freference));
%pval 0.0072
%corr 0.5450

%mk_dsrred
%can't subtract log values, so use these
nldsrAredt=metabolicgenesnoatp(konoatp==11180.2,:)';
nldsrBredt=metabolicgenesnoatp(konoatp==11181.2,:)';
dsrABredt=log(nldsrAredt + nldsrBredt);

dataX(18,:)=dsrABredt;
dataY(18,:)=log(h2s_freference);
%pval 0.8672
%corr 0.2258


dataX(19,:)=dsrABredt;
dataY(19,:)=log(srra_freference');
%pval 0.6432
%corr 0.2757


[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX',dataY');
writematrix(pval, "pval_rates_genes.txt");
writematrix(corr_obs, "corr_rates_genes.txt");

%Also, calculation the R2 values to compare to stepwise
mdlpsaB = fitlm(psaBt12, log(pprod_freference));
mdlHAO = fitlm(HAOt, log(nitri1_freference));
mdlamoAB = fitlm(amoABt, log(nitri1_freference));
mdlnosZ = fitlm(nosZt, log(dno2_freference+dno3_freference));
mdldsrABox = fitlm(dsrABoxt, log((soxo_freference + soxno2_freference + soxno3_freference)));
mdldsrABred = fitlm(dsrABredt, log(srra_freference'));
