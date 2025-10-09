clear("pval", "corr_obs", "crit_corr", "est_alpha", "seed_state", 'dataX', 'dataY');

%sig = 0
nlpetAt1=metabolicgenesnoatpHGnorm(konoatp==2634.1,:)';
nlpetAt2=metabolicgenesnoatpHGnorm(konoatp==2634.2,:)';
petAt12=log(nlpetAt1 + nlpetAt2);

dataX(1,:)=petAt12;
dataY(1,:)=log(metadata.CHLA);

%sig = 0
petBt=lnmetabolicgenesnoatp(konoatp==2635.1,:)';
nlpetBt1=metabolicgenesnoatpHGnorm(konoatp==2635.1,:)';
nlpetBt2=metabolicgenesnoatpHGnorm(konoatp==2635.2,:)';
petBt12=log(nlpetBt1 + nlpetBt2);

dataX(2,:)=petBt12;
dataY(2,:)=log(metadata.CHLA);

%NS
petCt=lnmetabolicgenesnoatp(konoatp==2636.1,:)';
nlpetCt1=metabolicgenesnoatpHGnorm(konoatp==2636.1,:)';
nlpetCt2=metabolicgenesnoatpHGnorm(konoatp==2636.2,:)';
petCt12=log(nlpetCt1 + nlpetCt2);

dataX(3,:)=petCt12;
dataY(3,:)=log(metadata.CHLA);

%sig=0
petDt=lnmetabolicgenesnoatp(konoatp==2637.1,:)';
nlpetDt1=metabolicgenesnoatpHGnorm(konoatp==2637.1,:)';
nlpetDt2=metabolicgenesnoatpHGnorm(konoatp==2637.2,:)';
petDt12=log(nlpetDt1 + nlpetDt2);

dataX(4,:)=petDt12;
dataY(4,:)=log(metadata.CHLA);

%NS
petEt=lnmetabolicgenesnoatp(konoatp==2638.1,:)';
nlpetEt1=metabolicgenesnoatpHGnorm(konoatp==2638.1,:)';
nlpetEt2=metabolicgenesnoatpHGnorm(konoatp==2638.2,:)';
petEt12=log(nlpetEt1 + nlpetEt2);

dataX(5,:)=petEt12;
dataY(5,:)=log(metadata.CHLA);

%NS
petFt=lnmetabolicgenesnoatp(konoatp==2639.1,:)';
nlpetFt1=metabolicgenesnoatpHGnorm(konoatp==2639.1,:)';
nlpetFt2=metabolicgenesnoatpHGnorm(konoatp==2639.2,:)';
petFt12=log(nlpetFt1 + nlpetFt2);

dataX(6,:)=petFt12;
dataY(6,:)=log(metadata.CHLA);

%NS
petGt=lnmetabolicgenesnoatp(konoatp==2640.1,:)';
nlpetGt1=metabolicgenesnoatpHGnorm(konoatp==2640.1,:)';
nlpetGt2=metabolicgenesnoatpHGnorm(konoatp==2640.2,:)';
petGt12=log(nlpetGt1 + nlpetGt2);

dataX(7,:)=petGt12;
dataY(7,:)=log(metadata.CHLA);

%NS
petHt=lnmetabolicgenesnoatp(konoatp==2641.1,:)';
nlpetHt1=metabolicgenesnoatpHGnorm(konoatp==2641.1,:)';
nlpetHt2=metabolicgenesnoatpHGnorm(konoatp==2641.2,:)';
petHt12=log(nlpetHt1 + nlpetHt2);

dataX(8,:)=petHt12;
dataY(8,:)=log(metadata.CHLA);

%NS
petJt=lnmetabolicgenesnoatp(konoatp==8906.1,:)';
nlpetJt1=metabolicgenesnoatpHGnorm(konoatp==8906.1,:)';
nlpetJt2=metabolicgenesnoatpHGnorm(konoatp==8906.2,:)';
petJt12=log(nlpetJt1 + nlpetJt2);

dataX(9,:)=petJt12;
dataY(9,:)=log(metadata.CHLA);

%pval=0
psaAt=lnmetabolicgenesnoatp(konoatp==2689.1,:)';
nlpsaAt1=metabolicgenesnoatpHGnorm(konoatp==2689.1,:)';
nlpsaAt2=metabolicgenesnoatpHGnorm(konoatp==2689.2,:)';
psaAt12=log(nlpsaAt1 + nlpsaAt2);

dataX(10,:)=psaAt12;
dataY(10,:)=log(metadata.CHLA);

%pval=0
psaBt=lnmetabolicgenesnoatp(konoatp==2690.1,:)';
nlpsaBt1=metabolicgenesnoatpHGnorm(konoatp==2690.1,:)';
nlpsaBt2=metabolicgenesnoatpHGnorm(konoatp==2690.2,:)';
psaBt12=log(nlpsaBt1 + nlpsaBt2);

dataX(11,:)=psaAt12;
dataY(11,:)=log(metadata.CHLA);

%ns
psaCt=lnmetabolicgenesnoatp(konoatp==2691.1,:)';
nlpsaCt1=metabolicgenesnoatpHGnorm(konoatp==2691.1,:)';
nlpsaCt2=metabolicgenesnoatpHGnorm(konoatp==2691.2,:)';
psaCt12=log(nlpsaCt1 + nlpsaCt2);

dataX(12,:)=psaCt12;
dataY(12,:)=log(metadata.CHLA);

%p=0.0388
psaDt=lnmetabolicgenesnoatp(konoatp==2692.1,:)';
nlpsaDt1=metabolicgenesnoatpHGnorm(konoatp==2692.1,:)';
nlpsaDt2=metabolicgenesnoatpHGnorm(konoatp==2692.2,:)';
psaDt12=log(nlpsaDt1 + nlpsaDt2);

dataX(13,:)=psaDt12;
dataY(13,:)=log(metadata.CHLA);

%pval=0.1804
psaEt=lnmetabolicgenesnoatp(konoatp==2693.1,:)';
nlpsaEt1=metabolicgenesnoatpHGnorm(konoatp==2693.1,:)';
nlpsaEt2=metabolicgenesnoatpHGnorm(konoatp==2693.2,:)';
psaEt12=log(nlpsaEt1 + nlpsaEt2);

dataX(14,:)=psaEt12;
dataY(14,:)=log(metadata.CHLA);

%p=0
psaFt=lnmetabolicgenesnoatp(konoatp==2694.1,:)';
nlpsaFt1=metabolicgenesnoatpHGnorm(konoatp==2694.1,:)';
nlpsaFt2=metabolicgenesnoatpHGnorm(konoatp==2694.2,:)';
psaFt12=log(nlpsaFt1 + nlpsaFt2);

dataX(15,:)=psaFt12;
dataY(15,:)=log(metadata.CHLA);

%NS
psbAt=lnmetabolicgenesnoatp(konoatp==2703.1,:)';
nlpsbAt1=metabolicgenesnoatpHGnorm(konoatp==2703.1,:)';
nlpsbAt2=metabolicgenesnoatpHGnorm(konoatp==2703.2,:)';
psbAt12=log(nlpsbAt1 + nlpsbAt2);

dataX(16,:)=psbAt12;
dataY(16,:)=log(metadata.CHLA);

%p=0
psbBt=lnmetabolicgenesnoatp(konoatp==2704.1,:)';
nlpsbBt1=metabolicgenesnoatpHGnorm(konoatp==2704.1,:)';
nlpsbBt2=metabolicgenesnoatpHGnorm(konoatp==2704.2,:)';
psbBt12=log(nlpsbBt1 + nlpsbBt2);

dataX(17,:)=psbBt12;
dataY(17,:)=log(metadata.CHLA);

%p=0
psbCt=lnmetabolicgenesnoatp(konoatp==2705.1,:)';
nlpsbCt1=metabolicgenesnoatpHGnorm(konoatp==2705.1,:)';
nlpsbCt2=metabolicgenesnoatpHGnorm(konoatp==2705.2,:)';
psbCt12=log(nlpsbCt1 + nlpsbCt2);

dataX(18,:)=psbCt12;
dataY(18,:)=log(metadata.CHLA);

%p=0.26
psbDt=lnmetabolicgenesnoatp(konoatp==2706.1,:)';
nlpsbDt1=metabolicgenesnoatpHGnorm(konoatp==2706.1,:)';
nlpsbDt2=metabolicgenesnoatpHGnorm(konoatp==2706.2,:)';
psbDt12=log(nlpsbDt1 + nlpsbDt2);

dataX(19,:)=psbDt12;
dataY(19,:)=log(metadata.CHLA);

%p=0
psbEt=lnmetabolicgenesnoatp(konoatp==2707.1,:)';
nlpsbEt1=metabolicgenesnoatpHGnorm(konoatp==2707.1,:)';
nlpsbEt2=metabolicgenesnoatpHGnorm(konoatp==2707.2,:)';
psbEt12=log(nlpsbEt1 + nlpsbEt2);

dataX(20,:)=psbEt12;
dataY(20,:)=log(metadata.CHLA);

%NS
psbFt=lnmetabolicgenesnoatp(konoatp==2708.1,:)';
nlpsbFt1=metabolicgenesnoatpHGnorm(konoatp==2708.1,:)';
nlpsbFt2=metabolicgenesnoatpHGnorm(konoatp==2708.2,:)';
psbFt12=log(nlpsbFt1 + nlpsbFt2);

dataX(21,:)=psbFt12;
dataY(21,:)=log(metadata.CHLA);

%NS
psbOt=lnmetabolicgenesnoatp(konoatp==2716.1,:)';
nlpsbOt1=metabolicgenesnoatpHGnorm(konoatp==2716.1,:)';
nlpsbOt2=metabolicgenesnoatpHGnorm(konoatp==2716.2,:)';
psbOt12=log(nlpsbOt1 + nlpsbOt2);

dataX(22,:)=psbOt12;
dataY(22,:)=log(metadata.CHLA);

[pval_chla, corr_obs_chla, crit_corr_chla, est_alpha_chla, seed_state_chla]=mult_comp_perm_corr(dataX',dataY');
writematrix(pval_chla, "pval_rates_genes_PS_chla.txt");
writematrix(corr_obs_chla, "corr_rates_genes_PS_chla.txt");