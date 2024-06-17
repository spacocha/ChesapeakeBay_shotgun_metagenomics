%run process_TPM_chesapeake_SPP
%load MATLAB_metadata_2.txt

MATLAB_metadata_2=readtable("MATLAB_metadata_2.txt");

%Add some stuff to MATLAB_metadata
MATLAB_metadata_3=MATLAB_metadata_2;

%log transform everything but AOU (negative)
%don't log PH, it's already log of H+
%MATLAB_metadata_3.logTEMP=log(MATLAB_metadata_2.TEMP);
%MATLAB_metadata_3.logSALT=log(MATLAB_metadata_2.SALT);
%MATLAB_metadata_3.logDO=log(MATLAB_metadata_2.DO);
%MATLAB_metadata_3.logPO4F=log(MATLAB_metadata_2.PO4F);
%MATLAB_metadata_3.logNO2F=log(MATLAB_metadata_2.NO2F);
%Not sure what to do about NO3F, there are 0 and negative numbers
%MATLAB_metadata_3.logNH4F=log(MATLAB_metadata_2.NH4F);
%MATLAB_metadata_3.logDON=log(MATLAB_metadata_2.DON);
%MATLAB_metadata_3.logDOP=log(MATLAB_metadata_2.DOP);
%MATLAB_metadata_3.logPN=log(MATLAB_metadata_2.PN);
%MATLAB_metadata_3.logPP=log(MATLAB_metadata_2.PP);
%MATLAB_metadata_3.logCHLA=log(MATLAB_metadata_2.CHLA);
%MATLAB_metadata_3.logPHEO=log(MATLAB_metadata_2.PHEO);
writetable(MATLAB_metadata_3, "MATLAB_metadata_3.txt");

indexnum=1;

%Make the merged value from the non-log and then do log after
petAnl=metabolicgenesnoatp(konoatp==2634.1,:)';
petBnl=metabolicgenesnoatp(konoatp==2635.1,:)';
petCnl=metabolicgenesnoatp(konoatp==2636.1,:)';
petDnl=metabolicgenesnoatp(konoatp==2637.1,:)';
petEnl=metabolicgenesnoatp(konoatp==2638.1,:)';
petFnl=metabolicgenesnoatp(konoatp==2639.1,:)';
petGnl=metabolicgenesnoatp(konoatp==2640.1,:)';
petHnl=metabolicgenesnoatp(konoatp==2641.1,:)';

allpet=log(petAnl + petBnl + petCnl + petDnl + petEnl + petFnl + petGnl + petHnl);

MATLAB_metadata_3.RESPONSE=allpet;
resultsallpet=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);

writetable(resultsallpet.Coefficients, 'resultsallpet.txt','WriteRowNames',true)

%FYI it's the same as resultsHAO.predict
scatter(resultsallpet.predict,allpet);

overallres(:,indexnum)=resultsallpet.Rsquared.Adjusted;
overallresnames(:,indexnum)="allpet";

indexnum=indexnum+1;

psaAnl=metabolicgenesnoatp(konoatp==2689.1,:)';
psaBnl=metabolicgenesnoatp(konoatp==2690.1,:)';
psaCnl=metabolicgenesnoatp(konoatp==2691.1,:)';
psaDnl=metabolicgenesnoatp(konoatp==2692.1,:)';
psaEnl=metabolicgenesnoatp(konoatp==2693.1,:)';
psaFnl=metabolicgenesnoatp(konoatp==2694.1,:)';

allpsa=log(psaAnl + psaBnl + psaCnl + psaDnl + psaEnl + psaFnl);

MATLAB_metadata_3.RESPONSE=allpsa;
resultsallpsa=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);

writetable(resultsallpsa.Coefficients, 'resultsallpsa.txt','WriteRowNames',true)

%FYI it's the same as resultsHAO.predict
scatter(resultsallpsa.predict,allpsa);

overallres(:,indexnum)=resultsallpsa.Rsquared.Adjusted;
overallresnames(:,indexnum)="allpsa";

indexnum=indexnum+1;

psbAnl=metabolicgenesnoatp(konoatp==2703.1,:)';
psbBnl=metabolicgenesnoatp(konoatp==2704.1,:)';
psbCnl=metabolicgenesnoatp(konoatp==2705.1,:)';
psbDnl=metabolicgenesnoatp(konoatp==2706.1,:)';
psbEnl=metabolicgenesnoatp(konoatp==2707.1,:)';
psbFnl=metabolicgenesnoatp(konoatp==2708.1,:)';
psbOnl=metabolicgenesnoatp(konoatp==2716.1,:)';

allpsb=log(psbAnl + psbBnl + psbCnl + psbDnl + psbEnl + psbFnl +psbOnl); 

MATLAB_metadata_3.RESPONSE=allpsb;
resultsallpsb=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);

writetable(resultsallpsb.Coefficients, 'resultsallpsb.txt','WriteRowNames',true)

%FYI it's the same as resultsHAO.predict
scatter(resultsallpsb.predict,allpsb);

overallres(:,indexnum)=resultsallpsb.Rsquared.Adjusted;
overallresnames(:,indexnum)="allpsb";

indexnum=indexnum+1;

%all together
allphoto = log(psbAnl + psbBnl + psbCnl + psbDnl + psbEnl + psbFnl +psbOnl + ...
psaAnl + psaBnl + psaCnl + psaDnl + psaEnl + psaFnl +...   )
petAnl + petBnl + petCnl + petDnl + petEnl + petFnl + petGnl + petHnl);

writematrix(allphoto, "allphoto.txt");

MATLAB_metadata_3.RESPONSE=allphoto;
resultsallphoto=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);

writetable(resultsallphoto.Coefficients, 'resultsallphoto.txt','WriteRowNames',true)

%FYI it's the same as resultsHAO.predict
scatter(resultsallphoto.predict,allphoto);

overallres(:,indexnum)=resultsallphoto.Rsquared.Adjusted;
overallresnames(:,indexnum)="allphoto";

indexnum=indexnum+1;

%run the regression model for EOF1
%add the response variable to the end
MATLAB_metadata_3.RESPONSE=uat(:,1);

resultsEOF1=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);

writetable(resultsEOF1.Coefficients, 'resultsEOF1.txt','WriteRowNames',true)

%FYI it's the same as resultsHAO.predict
scatter(resultsEOF1.predict,uat(:,1));
EOF1predictcorr=corr(resultsEOF1.predict,uat(:,1));

overallres(:,indexnum)=resultsEOF1.Rsquared.Adjusted;
overallresnames(:,indexnum)="EOF1";

indexnum=indexnum+1;

%run the regression model for EOF2
MATLAB_metadata_3.RESPONSE=uat(:,2);
resultsEOF2=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsEOF2.Coefficients, 'resultsEOF2.txt','WriteRowNames',true)

scatter(resultsEOF2.predict,uat(:,2));
EOF2predictcorr=corr(resultsEOF2.predict,uat(:,2));

overallres(:,indexnum)=resultsEOF2.Rsquared.Adjusted;
overallresnames(:,indexnum)="EOF2";

indexnum=indexnum+1;

%run the regression model for cEOF1
MATLAB_metadata_3.RESPONSE=uct(:,1);
resultscEOF1=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultscEOF1.Coefficients, 'resultscEOF1.txt','WriteRowNames',true)

scatter(resultscEOF1.predict,uct(:,1));

overallres(:,indexnum)=resultscEOF1.Rsquared.Adjusted;
overallresnames(:,indexnum)="cEOF1";

indexnum=indexnum+1;

%run the regression model for cEOF2
MATLAB_metadata_3.RESPONSE=uct(:,2);
resultscEOF2=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultscEOF2.Coefficients, 'resultscEOF2.txt','WriteRowNames',true)

scatter(resultscEOF2.predict,uct(:,2));

overallres(:,indexnum)=resultscEOF2.Rsquared.Adjusted;
overallresnames(:,indexnum)="cEOF2";

indexnum=indexnum+1;

%run the regression model for HAO
MATLAB_metadata_3.RESPONSE=HAOt;
resultsHAO=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsHAO.Coefficients, 'resultsHAO.txt','WriteRowNames',true)

%Check equation
scatter(resultsHAO.predict,HAOt);

%record results
overallres(:,indexnum)=resultsHAO.Rsquared.Adjusted;
overallresnames(:,indexnum)="HAO";

%advance index
indexnum=indexnum+1;

%Do the same with amoA
MATLAB_metadata_3.RESPONSE=amoAt;
resultsamoA=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsamoA.Coefficients, 'resultsamoA.txt','WriteRowNames',true)

%Same as resultsamoA.predict
scatter(resultsamoA.predict,amoAt);

overallres(:,indexnum)=resultsamoA.Rsquared.Adjusted;
overallresnames(:,indexnum)="amoA";

indexnum=indexnum+1;

%Do the same with amoB
MATLAB_metadata_3.RESPONSE=amoBt;
resultsamoB=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsamoB.Coefficients, 'resultsamoB.txt','WriteRowNames',true)

%Same as results.predict
scatter(resultsamoB.predict,amoBt);

overallres(:,indexnum)=resultsamoB.Rsquared.Adjusted;
overallresnames(:,indexnum)="amoB";

indexnum=indexnum+1;

%and amoC
MATLAB_metadata_3.RESPONSE=amoCt;
resultsamoC=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsamoC.Coefficients, 'resultsamoC.txt','WriteRowNames',true)

%Same as results.predict
scatter(resultsamoC.predict,amoCt);

overallres(:,indexnum)=resultsamoC.Rsquared.Adjusted;
overallresnames(:,indexnum)="amoC";
indexnum=indexnum+1;

%and narI
MATLAB_metadata_3.RESPONSE=narIt;
resultsnarI=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsnarI.Coefficients, 'resultsnarI.txt','WriteRowNames',true)

scatter(resultsnarI.predict,narIt);

overallres(:,indexnum)=resultsnarI.Rsquared.Adjusted;
overallresnames(:,indexnum)="narI";

indexnum=indexnum+1;

%and nirS
MATLAB_metadata_3.RESPONSE=nirSt;
resultsnirS=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsnirS.Coefficients, 'resultsnirS.txt','WriteRowNames',true)

%Check that the equation is right
scatter(resultsnirS.predict,nirSt);

%save results
overallres(:,indexnum)=resultsnirS.Rsquared.Adjusted;
overallresnames(:,indexnum)="nirS";

%advance the index
indexnum=indexnum+1;

%napA
MATLAB_metadata_3.RESPONSE=napAt;
resultsnapA=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsnapA.Coefficients, 'resultsnapA.txt','WriteRowNames',true)

%Check the equation
scatter(resultsnapA.predict,napAt);

%record results
overallres(:,indexnum)=resultsnapA.Rsquared.Adjusted;
overallresnames(:,indexnum)="napA";

%advance index
indexnum=indexnum+1;

%norB
MATLAB_metadata_3.RESPONSE=norBt;
resultsnorB=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsnorB.Coefficients, 'resultsnorB.txt','WriteRowNames',true)

%check equation
scatter(resultsnorB.predict,norBt);

%record
overallres(:,indexnum)=resultsnorB.Rsquared.Adjusted;
overallresnames(:,indexnum)="norB";

%advance index
indexnum=indexnum+1;

%nosZ
MATLAB_metadata_3.RESPONSE=nosZt;
resultsnosZ=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsnosZ.Coefficients, 'resultsnosZ.txt','WriteRowNames',true)

%check equation
scatter(resultsnosZ.predict,nosZt);

%record
overallres(:,indexnum)=resultsnosZ.Rsquared.Adjusted;
overallresnames(:,indexnum)="nosZ";

%advance index
indexnum=indexnum+1;

%dsrA
MATLAB_metadata_3.RESPONSE=dsrAoxt;
resultsdsrAox=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsdsrAox.Coefficients, 'resultsdsrAox.txt','WriteRowNames',true)

%check equation
scatter(resultsdsrAox.predict,dsrAoxt);

%record results
overallres(:,indexnum)=resultsdsrAox.Rsquared.Adjusted;
overallresnames(:,indexnum)="dsrAox";

%advance index
indexnum=indexnum+1;

%dsrB
MATLAB_metadata_3.RESPONSE=dsrBoxt;
resultsdsrBox=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsdsrBox.Coefficients, 'resultsdsrBox.txt','WriteRowNames',true)

%check equation
scatter(resultsdsrBox.predict,dsrBoxt);

%record results
overallres(:,indexnum)=resultsdsrBox.Rsquared.Adjusted;
overallresnames(:,indexnum)="dsrBox";

%advance index
indexnum=indexnum+1;

%hdra2
MATLAB_metadata_3.RESPONSE=hdrA2t;
resultshdrA2=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultshdrA2.Coefficients, 'resultshdrA2.txt','WriteRowNames',true)

%check equation
scatter(resultshdrA2.predict,hdrA2t);

%record results
overallres(:,indexnum)=resultshdrA2.Rsquared.Adjusted;
overallresnames(:,indexnum)="hdrA2";

%advance index
indexnum=indexnum+1;

%psbE
MATLAB_metadata_3.RESPONSE=psbEt;
resultspsbE=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspsbE.Coefficients, 'resultspsbE.txt','WriteRowNames',true)

%check equation
scatter(resultspsbE.predict,psbEt);

%record
overallres(:,indexnum)=resultspsbE.Rsquared.Adjusted;
overallresnames(:,indexnum)="psbE";

%advance index
indexnum=indexnum+1;


%psbB
MATLAB_metadata_3.RESPONSE=psaBt;
resultspsaB=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspsaB.Coefficients, 'resultspsaB.txt','WriteRowNames',true)

%check equation
scatter(resultspsaB.predict,psaBt);

%record results
overallres(:,indexnum)=resultspsaB.Rsquared.Adjusted;
overallresnames(:,indexnum)="psaB";

%advance index
indexnum=indexnum+1;


%psaA
MATLAB_metadata_3.RESPONSE=psaAt;
resultspsaA=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspsaA.Coefficients, 'resultspsaA.txt','WriteRowNames',true)

%check equation
scatter(resultspsaA.predict,psaAt);

%record results
overallres(:,indexnum)=resultspsaA.Rsquared.Adjusted;
overallresnames(:,indexnum)="psaA";

%advance index
indexnum=indexnum+1;

%psaC
MATLAB_metadata_3.RESPONSE=psaCt;
resultspsaC=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspsaC.Coefficients, 'resultspsaC.txt','WriteRowNames',true)

%check equation
scatter(resultspsaC.predict,psaCt);

%record results
overallres(:,indexnum)=resultspsaC.Rsquared.Adjusted;
overallresnames(:,indexnum)="psaC";

%advance index
indexnum=indexnum+1;

%petA
MATLAB_metadata_3.RESPONSE=petAt;
resultspetA=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspetA.Coefficients, 'resultspetA.txt','WriteRowNames',true)

%check equation
scatter(resultspetA.predict,petAt);

%record results
overallres(:,indexnum)=resultspetA.Rsquared.Adjusted;
overallresnames(:,indexnum)="petA";

%advance index
indexnum=indexnum+1;

%petB
MATLAB_metadata_3.RESPONSE=petBt;
resultspetB=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspetB.Coefficients, 'resultspetB.txt','WriteRowNames',true)

%check equation
scatter(resultspetB.predict,petBt);

%record results
overallres(:,indexnum)=resultspetB.Rsquared.Adjusted;
overallresnames(:,indexnum)="petB";

%advance index
indexnum=indexnum+1;

%petD
MATLAB_metadata_3.RESPONSE=petDt;
resultspetD=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspetD.Coefficients, 'resultspetD.txt','WriteRowNames',true)

%check equation
scatter(resultspetD.predict,petDt);

%record results
overallres(:,indexnum)=resultspetD.Rsquared.Adjusted;
overallresnames(:,indexnum)="petD";

%advance index
indexnum=indexnum+1;

%psa
MATLAB_metadata_3.RESPONSE=petDt;
resultspetD=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultspetD.Coefficients, 'resultspetD.txt','WriteRowNames',true)

%check equation
scatter(resultspetD.predict,petDt);

%record results
overallres(:,indexnum)=resultspetD.Rsquared.Adjusted;
overallresnames(:,indexnum)="petD";

%advance index
indexnum=indexnum+1;

%nifK
MATLAB_metadata_3.RESPONSE=nifKt;
resultsnifK=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsnifK.Coefficients, 'resultsnifK.txt','WriteRowNames',true)

%check equation
scatter(resultsnifK.predict,nifKt);

%record results
overallres(:,indexnum)=resultsnifK.Rsquared.Adjusted;
overallresnames(:,indexnum)="nifK";

%advance index
indexnum=indexnum + 1;

%soxA
MATLAB_metadata_3.RESPONSE=soxAt;
resultssoxA=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultssoxA.Coefficients, 'resultssoxA.txt','WriteRowNames',true)

%check equation
scatter(resultssoxA.predict,soxAt);

%record results
overallres(:,indexnum)=resultssoxA.Rsquared.Adjusted;
overallresnames(:,indexnum)="soxA";

%advance index
indexnum=indexnum+1;

%soxX
MATLAB_metadata_3.RESPONSE=soxXt;
resultssoxX=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultssoxX.Coefficients, 'resultssoxX.txt','WriteRowNames',true)

%check equation
scatter(resultssoxX.predict,soxXt);

%record results
overallres(:,indexnum)=resultssoxX.Rsquared.Adjusted;
overallresnames(:,indexnum)="soxX";

%advance index
indexnum=indexnum+1;

%soxC
MATLAB_metadata_3.RESPONSE=soxCt;
resultssoxC=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultssoxC.Coefficients, 'resultssoxC.txt','WriteRowNames',true)

%check equation
scatter(resultssoxC.predict,soxCt);

%record results
overallres(:,indexnum)=resultssoxC.Rsquared.Adjusted;
overallresnames(:,indexnum)="soxC";

%advance index
indexnum=indexnum+1;

%random housekeeping genes
%1870
IARSt=lnhousekeepinggenes(kohousekeeping==1870,:)';

%soxC
MATLAB_metadata_3.RESPONSE=IARSt;
resultsIARS=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsIARS.Coefficients, 'resultsIARS.txt','WriteRowNames',true)

%check equation
scatter(resultsIARS.predict,IARSt);

%record results
overallres(:,indexnum)=resultsIARS.Rsquared.Adjusted;
overallresnames(:,indexnum)="IARS";

%advance index
indexnum=indexnum+1;

%2601
nusGt=lnhousekeepinggenes(kohousekeeping==2601,:)';

MATLAB_metadata_3.RESPONSE=nusGt;
resultsnusG=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsnusG.Coefficients, 'resultsnusG.txt','WriteRowNames',true);

%check equation
scatter(resultsnusG.predict,nusGt);

%record results
overallres(:,indexnum)=resultsnusG.Rsquared.Adjusted;
overallresnames(:,indexnum)="nusG";

%advance index
indexnum=indexnum+1;

%1889
FARSAt=lnhousekeepinggenes(kohousekeeping==1889,:)';

MATLAB_metadata_3.RESPONSE=FARSAt;
resultsFARSA=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsFARSA.Coefficients, 'resultsFARSA.txt','WriteRowNames',true)

%check equation
scatter(resultsFARSA.predict,FARSAt);

%record results
overallres(:,indexnum)=resultsFARSA.Rsquared.Adjusted;
overallresnames(:,indexnum)="FARSA";

%advance index
indexnum=indexnum+1;

%Ribosomal genes seem to be better correlated with something
%2867
RPL11t=lnhousekeepinggenes(kohousekeeping==2867,:)';

MATLAB_metadata_3.RESPONSE=RPL11t;
resultsRPL11=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsRPL11.Coefficients, 'resultsRPL11.txt','WriteRowNames',true)

%check equation
scatter(resultsRPL11.predict,RPL11t);

%record results
overallres(:,indexnum)=resultsRPL11.Rsquared.Adjusted;
overallresnames(:,indexnum)="RPL11";

%advance index
indexnum=indexnum+1;

%2881
RPL18t=lnhousekeepinggenes(kohousekeeping==2881,:)';

MATLAB_metadata_3.RESPONSE=RPL18t;
resultsRPL18=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
writetable(resultsRPL18.Coefficients, 'resultsRPL18.txt','WriteRowNames',true)

%check equation
scatter(resultsRPL18.predict,RPL18t);

%recrd results
overallres(:,indexnum)=resultsRPL18.Rsquared.Adjusted;
overallresnames(:,indexnum)="RPL18";

%advance index
indexnum=indexnum+1;

writematrix(overallres,"overallres.txt");
writematrix(overallresnames,"overallresnames.txt");