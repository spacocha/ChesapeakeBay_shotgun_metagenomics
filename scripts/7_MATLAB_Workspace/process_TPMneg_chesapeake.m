load TPM_final_negrem_1.mat
load analyze_genes_2
ches_freference %Load comparison data from model
indenv=[1     2     3     4     5     6     7     8    10    11    12    13    15    16    17    22    23    24    25    26   28    29    30    31    32    33    34    35    36    37    38    39    40    41    42      45    46    47];
kogenes=TPMfinaltable{:,1};
TPM_genes_2_negrem=TPMfinaltable{:,indenv+1};
nonzeron=sum(TPM_genes_2_negrem'>0)';
kocommonn=kogenes(nonzeron>30);
kocommonn=kocommonn(1:end-1);

commongenesn=TPM_genes_2_negrem(nonzeron>30,:); %Pick genes found at at least 80% of sites
commongenesn=commongenesn(1:end-1,:); %Get rid of unassigned genes
minval=min(commongenesn(commongenesn(:)>0));
commongenesn(commongenesn(:)==0)=minval% Fill with a value that is smaller than any of the observed values
lncommongenesn=log(commongenesn);
lncommongenesanomn=lncommongenesn-mean(lncommongenesn')'; %Demean each set of genes

[uctn sctn vctn]=svd(lncommongenesanomn','econ');

ko=[362 363 366 367 368 370 371 372 374 376 380 381 390 392 394 395 ...
    958 1601 1602 2108 2109 2110 2111 2112 2113 2114 2115 2305 2567 2568 ... 
     2586 2588 2591 2634 2634 2635 2635 2636 2636 2637 2637 2638 2639 2640 2641 ...
    2689 2690 2691 2692 2693 2694 2703 2704 2705 2706 2707 2708 2716 3385 3388 3389 3390 ...
    4561 5301  8906 8928 8929 10534 10535 10944 10944.1 10944.2 10945.1 10945.2 10946 10946.1 10946.2 11180.1 11180.2 11181.1 11181.2 ...
    15864 15876 16257 16259 17222 17223 17224 17225 17226 17227 17229 17230 22622 23995];

konoatp=[362 363 366 367 368 370 371 372 374 376 380 381 390 392 394 395 ...
    958 1601 1602 2305 2567 2568 ... 
     2586 2588 2591 2634 2635 2636 2637 2638 2639 2640 2641 ...
    2689 2690 2691 2692 2693 2694 2703 2704 2705 2706 2707 2708 2716 3385 3388 3389 3390 ...
    4561 5301  8906 8928 8929 10534 10535 10944 10944.1 10944.2 10945.1 10945.2 10946.1 10946.2 11180.1 11180.2 11181.1 11181.2 ...
    15864 15876 16257 16259 17222 17223 17224 17225 17226 17227 17229 17230 22622 23995];

kohousekeeping=[1870 1872 1874 1876 1881 1887 1889 1890 1892 1937 2528 2601 2863 2864 2867 2871 2874 2876 2881 2886 2890 2906];


for j=1:length(ko)
    dummygenes=(TPM_genes_2_negrem(floor(kogenes)==ko(j),:));
    dummygenes2=TPM_genes_2_negrem(kogenes==ko(j),:);
    if (min(size(dummygenes))>2&isempty(dummygenes2))
      metabolicgenesn(j,:)=sum(dummygenes);
    elseif (min(size(dummygenes))==1&isempty(dummygenes2))
      metabolicgenesn(j,:)=dummygenes;
    end
    if ~isempty(dummygenes2)
    metabolicgenesn(j,:)=dummygenes2;
    end
end
for j=1:length(konoatp)
     if (sum(floor(kogenes)==konoatp(j))>1) 
         metabolicgenesnoatpn(j,:)=sum(TPM_genes_2_negrem(floor(kogenes)==konoatp(j),:));
     end
     if (sum(floor(kogenes)==konoatp(j))==1)
         metabolicgenesnoatpn(j,:)=TPM_genes_2_negrem(floor(kogenes)==konoatp(j),:);
     end
end
for j=1:length(kohousekeeping)
    if (sum(floor(kogenes)==kohousekeeping(j))>1)
      housekeepinggenesn(j,:)=sum(TPM_genes_2_negrem(floor(kogenes)==kohousekeeping(j),:));
    end
     if (sum(floor(kogenes)==kohousekeeping(j))==1)
      housekeepinggenesn(j,:)=(TPM_genes_2_negrem(floor(kogenes)==kohousekeeping(j),:));
    end
end
   meanhousekeepingn=mean(housekeepinggenesn);
    metabolicgenesn(metabolicgenesn(:)==0)=0.01;
    metabolicgenesnoatpn(metabolicgenesnoatpn(:)==0)=0.01;

    lnmetabolicgenesn=log(metabolicgenesn);
    lnmetabolicgenesanomn=lnmetabolicgenesn-mean(lnmetabolicgenesn')';
   

    lnmetabolicgenesnoatpn=log(metabolicgenesnoatpn);
    lnmetabolicgenesnoatpanomn=lnmetabolicgenesnoatpn-mean(lnmetabolicgenesnoatpn')';
    lnmetabolicgenesnoatpstdn=lnmetabolicgenesnoatpn-log(meanhousekeepingn);
    lnmetabolicgenesnoatpstdanomn=lnmetabolicgenesnoatpstdn-mean(lnmetabolicgenesnoatpstdn')';
    
[umtn smtn vmtn]=svd(lnmetabolicgenesanomn','econ');
[uatn sa_tn vatn]=svd(lnmetabolicgenesnoatpanomn','econ');
[ustdtn sstdtn vstdtn]=svd(lnmetabolicgenesnoatpstdanomn','econ');

[v1tn ind1tn]=sort(abs(corr(uatn(:,1),lnmetabolicgenesnoatpn')));
[v2tn ind2tn]=sort(abs(corr(uatn(:,2),lnmetabolicgenesnoatpn')));
[v3tn ind3tn]=sort(abs(corr(uatn(:,3),lnmetabolicgenesnoatpn')));
[v4tn ind4tn]=sort(abs(corr(uatn(:,4),lnmetabolicgenesnoatpn')));

[v1stdtn ind1stdtn]=sort(abs(corr(ustdtn(:,1),lnmetabolicgenesnoatpstdn')));
[v2stdtn ind2stdtn]=sort(abs(corr(ustdtn(:,2),lnmetabolicgenesnoatpstdn')));
[v3stdtn ind3stdtn]=sort(abs(corr(ustdtn(:,3),lnmetabolicgenesnoatpstdn')));
[v4stdtn ind4stdtn]=sort(abs(corr(ustdtn(:,4),lnmetabolicgenesnoatpstdn')));



nirBtn=lnmetabolicgenesn(ko==362,:)';
nirDtn=lnmetabolicgenesn(ko==363,:)';
nirAtn=lnmetabolicgenesn(ko==366,:)';
narBtn=lnmetabolicgenesn(ko==367,:)';
nirKtn=lnmetabolicgenesn(ko==368,:)';
narGZAtn=lnmetabolicgenesn(ko==370,:)';
narHYBtn=lnmetabolicgenesn(ko==371,:)';
nasAtn=lnmetabolicgenesn(ko==372,:)';
narIVtn=lnmetabolicgenesn(ko==374,:)';
nosZtn=lnmetabolicgenesn(ko==376,:)';
cysJtn=lnmetabolicgenesn(ko==380,:)';
cysItn=lnmetabolicgenesn(ko==381,:)';
cysHtn=lnmetabolicgenesn(ko==390,:)';
sirtn=lnmetabolicgenesn(ko==392,:)';
aprAtn=lnmetabolicgenesn(ko==394,:)';
aprBtn=lnmetabolicgenesn(ko==395,:)';
sattn=lnmetabolicgenesn(ko==958,:)';
rbcLtn=lnmetabolicgenesn(ko==1601,:)';
rbcStn=lnmetabolicgenesn(ko==1602,:)';
atpBtn=lnmetabolicgenesn(ko==2108,:)';
atpFtn=lnmetabolicgenesn(ko==2109,:)';
atpEtn=lnmetabolicgenesn(ko==2110,:)';
atpAtn=lnmetabolicgenesn(ko==2111,:)';
atpDtn=lnmetabolicgenesn(ko==2112,:)';
atpHtn=lnmetabolicgenesn(ko==2113,:)';
atpCtn=lnmetabolicgenesn(ko==2114,:)';
atpGtn=lnmetabolicgenesn(ko==2115,:)';
norCtn=lnmetabolicgenesn(ko==2305,:)';
napAtn=lnmetabolicgenesn(ko==2567,:)';
napBtn=lnmetabolicgenesn(ko==2568,:)';
nifDtn=lnmetabolicgenesn(ko==2586,:)';
nifHtn=lnmetabolicgenesn(ko==2588,:)';
nifKtn=lnmetabolicgenesn(ko==2591,:)';
petAtn=lnmetabolicgenesn(ko==2634,:)';
petBtn=lnmetabolicgenesn(ko==2635,:)';
petCtn=lnmetabolicgenesn(ko==2636,:)';
petDtn=lnmetabolicgenesn(ko==2637,:)';
petEtn=lnmetabolicgenesn(ko==2638,:)';
petFtn=lnmetabolicgenesn(ko==2639,:)';
petGtn=lnmetabolicgenesn(ko==2640,:)';
petHtn=lnmetabolicgenesn(ko==2641,:)';
psaAtn=lnmetabolicgenesn(ko==2689,:)';
psaBtn=lnmetabolicgenesn(ko==2690,:)';
psaCtn=lnmetabolicgenesn(ko==2691,:)';
psaDtn=lnmetabolicgenesn(ko==2692,:)';
psaEtn=lnmetabolicgenesn(ko==2693,:)';
psaFtn=lnmetabolicgenesn(ko==2694,:)';
psbAtn=lnmetabolicgenesn(ko==2703,:)';
psbBtn=lnmetabolicgenesn(ko==2704,:)';
psbCtn=lnmetabolicgenesn(ko==2705,:)';
psbDtn=lnmetabolicgenesn(ko==2706,:)';
psbEtn=lnmetabolicgenesn(ko==2707,:)';
psbFtn=lnmetabolicgenesn(ko==2708,:)';
psbOtn=lnmetabolicgenesn(ko==2716,:)';
nrfAtn=lnmetabolicgenesn(ko==3385,:)';
hdrA2tn=lnmetabolicgenesn(ko==3388,:)';
hdrB2tn=lnmetabolicgenesn(ko==3389,:)';
hdrC2tn=lnmetabolicgenesn(ko==3390,:)';
norBtn=lnmetabolicgenesn(ko==4561,:)';
sorAtn=lnmetabolicgenesn(ko==5301,:)';
petJtn=lnmetabolicgenesn(ko==8906,:)';
pufLtn=lnmetabolicgenesn(ko==8928,:)';
pufMtn=lnmetabolicgenesn(ko==8929,:)';
NRtn=lnmetabolicgenesn(ko==10534,:)';
HAOtn=lnmetabolicgenesn(ko==10535,:)';
pmoAamoAtn=lnmetabolicgenesn(ko==10944,:)';
amoAtn=lnmetabolicgenesn(ko==10944.1,:)';
pmoAtn=lnmetabolicgenesn(ko==10944.2,:)';
pmoBamoBtn=lnmetabolicgenesn(ko==10945,:)';
amoBtn=lnmetabolicgenesn(ko==10945.1,:)';
pmoBtn=lnmetabolicgenesn(ko==10945.2,:)';
pmoCamoCtn=lnmetabolicgenesn(ko==10946,:)';
amoCtn=lnmetabolicgenesn(ko==10946.1,:)';
pmoCtn=lnmetabolicgenesn(ko==10946.2,:)';
dsrAtnO=lnmetabolicgenesn(ko==11180.1,:)';
dsrAtnR=lnmetabolicgenesn(ko==11180.2,:)';
dsrBtnO=lnmetabolicgenesn(ko==11181.1,:)';
dsrBtnR=lnmetabolicgenesn(ko==11181.2,:)';
nirStn=lnmetabolicgenesn(ko==15864,:)';
nrfHtn=lnmetabolicgenesn(ko==15876,:)';
mxaCtn=lnmetabolicgenesn(ko==16257,:)';
mxaLtn=lnmetabolicgenesn(ko==16259,:)';
soxAtn=lnmetabolicgenesn(ko==17222,:)';
soxXtn=lnmetabolicgenesn(ko==17223,:)';
soxBtn=lnmetabolicgenesn(ko==17224,:)';
soxCtn=lnmetabolicgenesn(ko==17225,:)';
soxYtn=lnmetabolicgenesn(ko==17226,:)';
soxZtn=lnmetabolicgenesn(ko==17227,:)';
fccBtn=lnmetabolicgenesn(ko==17229,:)';
fccAtn=lnmetabolicgenesn(ko==17230,:)';
soxDtn=lnmetabolicgenesn(ko==22622,:)';
xoxFtn=lnmetabolicgenesn(ko==23995,:)';

names=['  nirB  '
       '  nirD  '
       '  nirA  '
       '  narB  '
       '  nirK  '
       ' narGZA '
       ' narHYB '
       '  nasA  '
       ' narIV  '
       '  nosZ  '
       '  cysJ  '
       '  cysI  '
       '  cysH  '
       '  sir   '
       '  aprA  '
       '  aprB  '
       '  sat   '
       '  rbcL  '
       '  rbcS  '
  %     '  atpB  '
  %     '  atpF  '
  %     '  atpE  '
  %     '  atpA  '
  %     '  atpD  '
  %     '  atpH  '
  %     '  atpC  '
  %     '  atpG  '
       '  norC  '
       '  napA  '
       '  napB  '
       '  nifD  '
       '  nifH  '
       '  nifK  '
       '  petA  '
       '  petB  '
       '  petC  '
       '  petD  '
       '  petE  '
       '  petF  '
       '  petG  '
       '  petH  '
       '  psaA  '
       '  psaB  '
       '  psaC  '
       '  psaD  '
       '  psaE  '
       '  psaF  '
       '  psbA  '
       '  psbB  '
       '  psbC  '
       '  psbD  '
       '  psbE  '
       '  psbF  '
       '  psbO  '
       '  nrfA  '
       '  hdrA2 '
       '  hdrB2 '
       '  hdrC2 '
       '  norB  '
       '  sorA  '
       '  petJ  '
       '  pufL  '
       '  pufM  '
       '    NR  '
       '   HAO  '
       'pmoAamoA'
       '  amoA  '
       '  pmoA  '
       'pmoBamoB'
       '  amoB  '
       '  pmoB  '
       'pmoCamoC'
       '  amoC  '
       '  pmoC  '
       '  dsrAO '
       '  drsAR '
       '  dsrBO '
       '  dsrBR '
       '  nirS  '
       '  nrfH  '
       '  mxaC  '
       '  mxaL  '
       '  soxA  '
       '  soxX  '
       '  soxB  '
       '  soxC  '
       '  soxY  '
       '  soxZ  '
       '  fccB  '
       '  fccA  '
       '  soxD  '
       '  xoxF  '];
   
close all
%% 

figure(1)
subplot(321)
scatter(log(nh4_freference),log(nitri1_freference),40,Yearday,'filled');caxis([100 250])
title('(a) Modeled Q_{nitri1}, vs. NH_4, Time')
xlabel('ln NH_4 in mmol m^{-3}')
ylabel('ln Rate in mmol m^{-3} d^{-1}')


subplot(322)
scatter(log(NH4F),HAOtn,40,Yearday,'filled');caxis([100 250]);
title('(b) Observed HAO vs. NH_4, Time')
ylabel('ln HAO in TPM')
xlabel('ln NH_4')

subplot(323)
scatter(HAOtn,log(nh4_freference),40,Yearday,'filled');caxis([100 250])
title('(c) Modelled Q_{nitri1}, vs. Obs.HAO,Time')
xlabel('ln HAO in TPM')
ylabel('ln Rate in mol m^{-3} d^{-1}')

subplot(324)
scatter(log(exp(amoAtn)+exp(pmoAamoAtn)-minval),log(nh4_freference),40,Yearday,'filled');caxis([100 250])
title('(d) Modelled Q_{nitri1}, vs. Obs.amoA,Time')
xlabel('ln amoA in TPM')
ylabel('ln Rate in mol m^{-3} d^{-1}')

subplot(325)
scatter(amoBtn,log(nh4_freference),40,Yearday,'filled');caxis([100 250])
title('(e) Modelled Q_{nitri1} vs. Obs.amoB, Time')
xlabel('ln amoB in TPM')
ylabel('ln Rate mol m^{-3} d^{-1}')
text(-2,-5,['Colors:Yearday'])
subplot(326)
scatter(log(exp(amoCtn)+exp(pmoCamoCtn)-minval),log(nh4_freference),40,Yearday,'filled');hc=colorbar;caxis([100 250]);set(hc,'Location','south')
title('(f) Modelled Q_{nitri1}, vs. Obs. amoC, Time')
xlabel('ln amoC in TPM')
ylabel('ln Rate in mol m^{-3} d^{-1}')
%% 
figure(5)
subplot(231)
scatter(log(CHLA),log(sum(TPMfinaltable{floor(TPMfinaltable{:,1})==2634,indenv+1})'+.01),40, Yearday,'filled')
title('(a) Obs. Chl vs. {\it petA}')
xlabel('log Chl (mg/m^3)')
ylabel('log TPM')
grid on
axis('square')

subplot(232)
scatter(log(chl_freference),log(sum(TPMfinaltable{floor(TPMfinaltable{:,1})==2634,indenv+1})'+.01),40, Yearday,'filled')
title('(b) Mod. Chl vs. {\it petA}')
xlabel('log Chl (mg/m^3)')
ylabel('log TPM')
grid on
axis('square')

subplot(233)
scatter(log(pprod_freference),log(sum(TPMfinaltable{floor(TPMfinaltable{:,1})==2634,indenv+1})'+.01),40, Yearday,'filled')
title('(c) Mod. Prim. Prod vs. {\it petA}')
xlabel('log Prim. Prod. (mmol N/m^3/day)')
ylabel('log TPM')
grid on
axis('square')

subplot(234)
scatter(log(CHLA),log((TPMfinaltable{(TPMfinaltable{:,1})==2634.1,indenv+1})'+.01),40, Yearday,'filled')
title('(c) Obs. Chl vs. Cyano {\it petA}')
xlabel('log Chl (mg/m^3)')
ylabel('log TPM')
grid on
axis('square')

subplot(235)
scatter(log(chl_freference),log((TPMfinaltable{(TPMfinaltable{:,1})==2634.1,indenv+1})'+.01),40, Yearday,'filled')
title('(e) Mod. Chl vs. Cyano{\it petA}')
xlabel('log Chl (mg/m^3)')
ylabel('log TPM')
grid on
hl=colorbar('Location','South');caxis([100 250])
axis('square')
subplot(236)
scatter(log(pprod_freference),log((TPMfinaltable{(TPMfinaltable{:,1})==2634.1,indenv+1})'+.01),40, Yearday,'filled')
title('(c) Mod. Prim. Prod vs. Cyano {\it petA}')
xlabel('log Prim. Prod. (mmol N/m^3/day)')
ylabel('log TPM')
grid on
axis('square')
%% 
figure(5)
subplot(231)
scatter(log(CHLA),log(sum(TPMfinaltable{floor(TPMfinaltable{:,1})==2690,indenv+1})'+.01),40, Yearday,'filled')
title('(a) Obs. Chl vs. {\it psaB}')
xlabel('log Chl (mg/m^3)')
ylabel('log nRPMK')
grid on
axis('square')

subplot(232)
scatter(log(chl_freference),log(sum(TPMfinaltable{floor(TPMfinaltable{:,1})==2690,indenv+1})'+.01),40, Yearday,'filled')
title('(b) Mod. Chl vs. {\it psaB}')
xlabel('log Chl (mg/m^3)')
ylabel('log nRPMK')
grid on
axis('square')

subplot(233)
scatter(log(pprod_freference),log(sum(TPMfinaltable{floor(TPMfinaltable{:,1})==2690,indenv+1})'+.01),40, Yearday,'filled')
title('(c) Mod. Prim. Prod vs. {\it psaB}')
xlabel('log Prim. Prod. (mmol N/m^3/day)')
ylabel('log nRPMK')
grid on
axis('square')

subplot(234)
scatter(log(CHLA),log((TPMfinaltable{(TPMfinaltable{:,1})==2690.1,indenv+1})'+.01),40, Yearday,'filled')
title('(c) Obs. Chl vs. Cyano {\it psaB}')
xlabel('log Chl (mg/m^3)')
ylabel('log nRPMK')
grid on
axis('square')

subplot(235)
scatter(log(chl_freference),log((TPMfinaltable{(TPMfinaltable{:,1})==2690.1,indenv+1})'+.01),40, Yearday,'filled')
title('(e) Mod. Chl vs. Cyano{\it psaB}')
xlabel('log Chl (mg/m^3)')
ylabel('log nRPMK')
grid on
hl=colorbar('Location','South');caxis([100 250])
axis('square')

subplot(236)
scatter(log(pprod_freference),log((TPMfinaltable{(TPMfinaltable{:,1})==2690.1,indenv+1})'+.01),40, Yearday,'filled')
title('(f) Mod. Prim. Prod vs. Cyano {\it psaB}')
xlabel('log Prim. Prod. (mmol N/m^3/day)')
ylabel('log nRPMK')
grid on
axis('square')