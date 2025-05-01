%Can we predict any of the genes with the model variables
%run process_TPM_chesapeake_SPP and ches_freference_SPP

nlpsaBt1=metabolicgenesnoatp(konoatp==2690.1,:)';
nlpsaBt2=metabolicgenesnoatp(konoatp==2690.2,:)';
psaBt12=log(nlpsaBt1 + nlpsaBt2);

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
%Missing PP
MATLAB_modeldata_3.CHL=chl_freference;
%Missing PHEO
MATLAB_modeldata_3.pH=pH_freference';

MATLAB_modeldata_4 = struct2table(MATLAB_modeldata_3);

MATLAB_modeldata_4.RESPONSE=psaBt12;
resultspsaB=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);

MATLAB_modeldata_4.RESPONSE=HAOt;
resultsHAO=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);

MATLAB_modeldata_4.RESPONSE=amoABt;
resultsamoA=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);

MATLAB_modeldata_4.RESPONSE=nosZt;
resultsnosZ=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);



MATLAB_modeldata_4.RESPONSE=dsrABoxt;
resultsdsrABox=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);

MATLAB_modeldata_4.RESPONSE=dsrABredt;
resultsdsrABred=stepwiselm(MATLAB_modeldata_4,'PEnter',0.01);
