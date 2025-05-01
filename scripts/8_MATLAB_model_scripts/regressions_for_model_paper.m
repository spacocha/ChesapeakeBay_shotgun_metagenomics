%run process_TPM_chesapeake_SPP
%Run ches_freference_SPP.m to get the modeled data

%mk_psaB
%can't subtract log values, so use these
nlpsaBt1=metabolicgenesnoatp(konoatp==2690.1,:)';
nlpsaBt2=metabolicgenesnoatp(konoatp==2690.2,:)';
psaBt12=log(nlpsaBt1 + nlpsaBt2);

%make amoAB
nlamoAt=metabolicgenesnoatp(konoatp==10944.1,:)';
nlamoBt=metabolicgenesnoatp(konoatp==10945.1,:)';
amoABt=log(nlamoAt + nlamoBt);

%mk_dsrABox
%can't subtract log values, so use these
nldsrAoxt=metabolicgenesnoatp(konoatp==11180.1,:)';
nldsrBoxt=metabolicgenesnoatp(konoatp==11181.1,:)';
dsrABoxt=log(nldsrAoxt + nldsrBoxt);

%mk_dsrABred
%can't subtract log values, so use these
nldsrAredt=metabolicgenesnoatp(konoatp==11180.2,:)';
nldsrBredt=metabolicgenesnoatp(konoatp==11181.2,:)';
dsrABredt=log(nldsrAredt + nldsrBredt);

%Calculate the adjusted R2 values for model
mdlpsaB = fitlm(psaBt12, log(pprod_freference));
results(1,1)=mdlpsaB.Rsquared.Adjusted;
mdlHAO = fitlm(HAOt, log(nitri1_freference));
results(2,1)=mdlHAO.Rsquared.Adjusted;
mdlamoAB = fitlm(amoABt, log(nitri1_freference));
results(3,1)=mdlamoAB.Rsquared.Adjusted;
mdlnosZ = fitlm(nosZt, log(dno2_freference+dno3_freference));
results(4,1)=mdlnosZ.Rsquared.Adjusted;
mdldsrABox = fitlm(dsrABoxt, log((soxo_freference + soxno2_freference + soxno3_freference)));
results(5,1)=mdldsrABox.Rsquared.Adjusted;
mdldsrABred = fitlm(dsrABredt, log(srra_freference'));
results(6,1)=mdldsrABred.Rsquared.Adjusted;

%Calculate the adjusted R2 values for modeled variables
MATLAB_modeldata_3.pH=pH_freference';
MATLAB_modeldata_3.TEMP=temp_freference;
MATLAB_modeldata_3.SALT=salt_freference;
MATLAB_modeldata_3.DO=DO_freference;
MATLAB_modeldata_3.PO4=po4_freference;
MATLAB_modeldata_3.NO2=no2_freference;
MATLAB_modeldata_3.NO3=no3_freference;
MATLAB_modeldata_3.NH4=nh4_freference;
MATLAB_modeldata_3.DON=semilabile_don_freference;
MATLAB_modeldata_3.DOP=semilabile_dop_freference;
MATLAB_modeldata_3.PN=refract_don_freference;
%Don't include pigments
%MATLAB_modeldata_3.CHL=chl_freference;
%Missing PHEO

MATLAB_modeldata_4 = struct2table(MATLAB_modeldata_3);

MATLAB_modeldata_4.RESPONSE=psaBt12;
resultspsaB=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);
results(1,2)=resultspsaB.Rsquared.Adjusted;

MATLAB_modeldata_4.RESPONSE=HAOt;
resultsHAO=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);
results(2,2)=resultsHAO.Rsquared.Adjusted;

MATLAB_modeldata_4.RESPONSE=amoABt;
resultsamoA=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);
results(3,2)=resultsamoA.Rsquared.Adjusted;

MATLAB_modeldata_4.RESPONSE=nosZt;
resultsnosZ=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);
results(4,2)=resultsnosZ.Rsquared.Adjusted;

MATLAB_modeldata_4.RESPONSE=dsrABoxt;
resultsdsrABox=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);
results(5,2)=resultsdsrABox.Rsquared.Adjusted;

MATLAB_modeldata_4.RESPONSE=dsrABredt;
resultsdsrABred=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);
results(6,2)=resultsdsrABred.Rsquared.Adjusted;

%Calculate adjusted R2 for observations
%remove AOU (calulcated) and CHLA (biological)
MATLAB_metadata_2=readtable("MATLAB_metadata_2.txt");
MATLAB_metadata_3=removevars(MATLAB_metadata_2,{'AOU'});
MATLAB_metadata_3=removevars(MATLAB_metadata_3,{'CHLA'});
MATLAB_metadata_3=removevars(MATLAB_metadata_3,{'PHEO'});
%Remove PP (not in model)
MATLAB_metadata_3=removevars(MATLAB_metadata_3,{'PP'});

MATLAB_metadata_3.RESPONSE=psaBt12;
resultspsaBt12=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
results(1,3)=resultspsaBt12.Rsquared.Adjusted;

MATLAB_metadata_3.RESPONSE=psaBt;
resultspsaB=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);

MATLAB_metadata_3.RESPONSE=HAOt;
resultsHAOt=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
results(2,3)=resultsHAOt.Rsquared.Adjusted;

MATLAB_metadata_3.RESPONSE=amoABt;
resultsamoABt=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
results(3,3)=resultsamoABt.Rsquared.Adjusted;

MATLAB_metadata_3.RESPONSE=nosZt;
resultsnosZt=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
results(4,3)=resultsnosZt.Rsquared.Adjusted;

MATLAB_metadata_3.RESPONSE=dsrABoxt;
resultsdsrABoxt=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
results(5,3)=resultsdsrABoxt.Rsquared.Adjusted;

MATLAB_metadata_3.RESPONSE=dsrABredt;
resultsdsrABredt=stepwiselm(MATLAB_metadata_3,'PEnter',0.01);
results(6,3)=resultsdsrABredt.Rsquared.Adjusted;

%Make the summary table related to the potential for rates
%to influence gene abundance
results(:,4)=results(:,1).*results(:,3)./results(:,2);
writematrix(results, "Rsquared_results.txt");