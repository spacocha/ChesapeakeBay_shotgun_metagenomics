%% 
%load TPM_final_table.KO.neg.merged.divided.final.txt
TPM_final_table_KO_neg_merged_divided_final=readtable("TPM_final_table.KO.neg.merged.divided.final.HAO.txt", 'TrimNonNumeric', true);

%limit to only valid headers
indenv=[1     2     3     4     5     6     7     8    10    11    12    13    15    16    17    22    23    24    25    26   28    29    30    31    32    33    34    35    36    37    38    39    40    41    42      45    46    47];
kogenes=TPM_final_table_KO_neg_merged_divided_final{:,1};
TPM_genes_2_NR=TPM_final_table_KO_neg_merged_divided_final{:,indenv+1};

%randomize
[m,n]=size(TPM_genes_2_NR);
[~,I]=sort(rand([m,n]));
J=repmat(1:n,m,1);
TPM_genes_2 = TPM_genes_2_NR(sub2ind([m,n],I,J));

nonzero=sum(TPM_genes_2'>0)';
%Choose genes that are non-zero across more than 30 (80%) of
%samples
kocommon=kogenes(nonzero>30);
kocommon=kocommon(1:end-1);

commongenes=TPM_genes_2(nonzero>30,:); %Pick genes found at at least 80% of sites
commongenes=commongenes(1:end-1,:); %Get rid of unassigned genes
minval=min(commongenes(commongenes(:)>0));
commongenes(commongenes(:)==0)=0.01; % Fill with a value that is smaller than any of the observed values but check that with minval
lncommongenes=log(commongenes);
lncommongenesanom=lncommongenes-mean(lncommongenes')'; %Demean each set of genes

[uct sct vct]=svd(lncommongenesanom','econ');

%determine the amount of variation for each EOF
comEOF1=sct(1,1).^2/sum(diag(sct.^2));
comEOF2=sct(2,2).^2/sum(diag(sct.^2));
comEOF3=sct(3,3).^2/sum(diag(sct.^2));
comEOF4=sct(4,4).^2/sum(diag(sct.^2));

%Check the correlation with salt after importing MATLAB_metadata.txt
%cEOF_SALT_corr=corr(MATLAB_metadata.SALT,uct(:,:));
%max(abs(cEOF_SALT_corr));

%
%ko=[362 363 366 367 368 370 371 372 374 376 380 381 390 392 394 395 ...
%    958 1601 1602 2108 2109 2110 2111 2112 2113 2114 2115 2305 2567 2568 ... 
%     2586 2588 2591 2634 2635 2636 2637 2638 2639 2640 2641 ...
%    2689 2690 2691 2692 2693 2694 2703 2704 2705 2706 2707 2708 2716 3385 3388 3389 3390 ...
%    4561 5301  8906 8928 8929 10534 10535 10944 10944.1 10944.2 10945.1 10945.2 10946.1 10946.2 11180.1 11180.2 11181.1 11181.2 ...
%    15864 15876 16257 16259 17222 17223 17224 17225 17226 17227 17229 17230 22622 23995];

konoatp=[362 363 366 367 368 370 371 372 374 376 380 381 390 392 394 395 ...
    958 1601 1602 2305 2567 2568 ... 
    2586 2588 2591 2634.1 2634.2 2635.1 2635.2 2636.1 2636.2 2637.1 2637.2 2638.1 2638.2 2639.1 2639.2 2640.1 2640.2 2641.1 2641.2 ...
    2689.1 2689.2 2690.1 2690.2 2691.1 2691.2 2692.1 2692.2 2693.1 2693.2 2694.1 2694.2 ...
    2703.1 2703.2 2704.1 2704.2 2705.1 2705.2 2706.1 2706.2 2707.1 2707.2 2708.1 2708.2 2716.1 2716.2 3385 3388 3389 3390 ...
    4561 5301  8906.1 8906.2 8928 8929 10534 10535.1 10535.2 10944.1 10944.2 10945.1 10945.2 10946.1 10946.2 11180.1 11180.2 11181.1 11181.2  ...
    15864 15876 16257 16259 17222 17223 17224 17225 17226 17227 17229 17230 22622 23995];


%single copy genes from non-model paper
%kohousekeeping=[1870 1872 1874 1876 1881 1887 1889 1890 1892 1937 2528 2601 2863 2864 2867 2871 2874 2876 2881 2886 2890 2906];

%A new set, not sure where it came frmo
%kohousekeeping=[19788 1889 1887 1875 1883 1869 1873 1409 20541 3110];

%single copy genes from Dupont et al The ISME Journal (2012) 6, 1186–1199
%Table S1, used in normalizing gene coverage as in https://doi.org/10.1128/AEM.02168-17.
%all yield essentially the same results
kohousekeeping =[2950 2890 2988 604 859 927 942 1866 1868 1869 1870 1872 1873 1875 1880 1881 1883 1887 1889 1890 1892 1937 2469 2470 ...
    2835 2838 2863 2864 2867 2871 2876 2892 2895 2906 2956 2961 2996 3106 3110 3595 3685 3687 4043 6942 14164 2992 3043 3545 2874 2967 ...
    2519 2520 2601 3075 3664 2948 2933 2886 2878 3977 2952 2904 2316 3702 9759 2879 2884 2887 2888 2899 2902 2914 2916 2926 2935 2939 ...
    2959 3979 3040 3046 3438 3596 3070 4075 2994 2968 2834 3553 2982 2357 2990 2986 2313 2338 2343 7042 3076 2911 2963 566 1972 2600 ...
    2965 2881 2946];

%for j=1:length(ko)
%    metabolicgenes(j,:)=TPM_genes_2(kogenes==ko(j),:);
%end
for j=1:length(konoatp)
     metabolicgenesnoatp(j,:)=TPM_genes_2(kogenes==konoatp(j),:);
end
for j=1:length(kohousekeeping)
      housekeepinggenes(j,:)=TPM_genes_2(kogenes==kohousekeeping(j),:);
end
   meanhousekeeping=mean(housekeepinggenes);
    %metabolicgenes(metabolicgenes(:)==0)=0.01;
    
    minvalnoatp=min(metabolicgenesnoatp(metabolicgenesnoatp(:)>0));
    %Set 0 values to 0.01, but check that it's lower than minvalnoatp;
    %add a pseudocount to 0 values
    metabolicgenesnoatp(metabolicgenesnoatp(:)==0)=0.01;

    %add a pseudocount to 0 values
    minvalhouse=min(housekeepinggenes(housekeepinggenes (:)>0));
    housekeepinggenes(housekeepinggenes(:)==0)=0.01;

    %lnmetabolicgenes=log(metabolicgenes);
    %lnmetabolicgenesanom=lnmetabolicgenes-mean(lnmetabolicgenes')';
    
    %normalize based on mean housekeeping gene abundance
    metabolicgenesnoatpHGnorm=metabolicgenesnoatp./meanhousekeeping;
    housekeepinggenesHGnorm=housekeepinggenes./meanhousekeeping;

    lnmetabolicgenesnoatp=log(metabolicgenesnoatpHGnorm);
    lnhousekeepinggenes=log(housekeepinggenesHGnorm);

    %demean the distributions for metabolic genes
    lnmetabolicgenesnoatpanom=lnmetabolicgenesnoatp-mean(lnmetabolicgenesnoatp')';

    %demean the distribution of housekeeping genes
    lnhousekeepinggenesnom=lnhousekeepinggenes-mean(lnhousekeepinggenes')';
    meanlnhousekeepingnom=mean(lnhousekeepinggenesnom);
    uctmlhncorr=corr(uct(:,1),meanlnhousekeepingnom');
    %This shows a good relationship if I want to make a figure
    %scatter(uct(:,1),meanlmhousekeepingnom')

    %not standardizing to meanhousekeeping genes anymore
    %lnmetabolicgenesnoatpstd=lnmetabolicgenesnoatp-log(meanhousekeeping);
    %lnmetabolicgenesnoatpstdanom=lnmetabolicgenesnoatpstd-mean(lnmetabolicgenesnoatpstd')';
    
%[umt smt vmt]=svd(lnmetabolicgenesanom','econ');
[uat sa_t vat]=svd(lnmetabolicgenesnoatpanom','econ');

%Check the correlation with salt after importing MATLAB_metadata.txt
%mEOF_SALT_corr=corr(MATLAB_metadata.SALT,uat(:,:));
%max(abs(mEOF_SALT_corr))

%percent of variance explained by each EOF
varEOF1=sa_t(1,1).^2/sum(diag(sa_t.^2));
varEOF2=sa_t(2,2).^2/sum(diag(sa_t.^2));
varEOF3=sa_t(3,3).^2/sum(diag(sa_t.^2));
varEOF4=sa_t(4,4).^2/sum(diag(sa_t.^2));
varEOF5=sa_t(5,5).^2/sum(diag(sa_t.^2));

%[ustdt sstdt vstdt]=svd(lnmetabolicgenesnoatpstdanom','econ');

%Not used. Instead, use mult_comp_perm_corr.m 
%described at the bottom
%[v1t ind1t]=sort(abs(corr(uat(:,1),lnmetabolicgenesnoatp')));
%[v2t ind2t]=sort(abs(corr(uat(:,2),lnmetabolicgenesnoatp')));
%[v3t ind3t]=sort(abs(corr(uat(:,3),lnmetabolicgenesnoatp')));
%[v4t ind4t]=sort(abs(corr(uat(:,4),lnmetabolicgenesnoatp')));

%[v1stdt ind1stdt]=sort(abs(corr(ustdt(:,1),lnmetabolicgenesnoatpstd')));
%[v2stdt ind2stdt]=sort(abs(corr(ustdt(:,2),lnmetabolicgenesnoatpstd')));
%[v3stdt ind3stdt]=sort(abs(corr(ustdt(:,3),lnmetabolicgenesnoatpstd')));
%[v4stdt ind4stdt]=sort(abs(corr(ustdt(:,4),lnmetabolicgenesnoatpstd')));



nirBt=lnmetabolicgenesnoatp(konoatp==362,:)';
nirDt=lnmetabolicgenesnoatp(konoatp==363,:)';
nirAt=lnmetabolicgenesnoatp(konoatp==366,:)';
narBt=lnmetabolicgenesnoatp(konoatp==367,:)';
nirKt=lnmetabolicgenesnoatp(konoatp==368,:)';
nirSt=lnmetabolicgenesnoatp(konoatp==15864,:)';
narGt=lnmetabolicgenesnoatp(konoatp==370,:)';
narHt=lnmetabolicgenesnoatp(konoatp==371,:)';
nasCt=lnmetabolicgenesnoatp(konoatp==372,:)';
narIt=lnmetabolicgenesnoatp(konoatp==374,:)';
nosZt=lnmetabolicgenesnoatp(konoatp==376,:)';
cysJt=lnmetabolicgenesnoatp(konoatp==380,:)';
cysIt=lnmetabolicgenesnoatp(konoatp==381,:)';
cysHt=lnmetabolicgenesnoatp(konoatp==390,:)';
sirt=lnmetabolicgenesnoatp(konoatp==392,:)';
aprAt=lnmetabolicgenesnoatp(konoatp==394,:)';
aprBt=lnmetabolicgenesnoatp(konoatp==395,:)';
satt=lnmetabolicgenesnoatp(konoatp==958,:)';
rbcLt=lnmetabolicgenesnoatp(konoatp==1601,:)';
rbcSt=lnmetabolicgenesnoatp(konoatp==1602,:)';
norCt=lnmetabolicgenesnoatp(konoatp==2305,:)';
napAt=lnmetabolicgenesnoatp(konoatp==2567,:)';
napBt=lnmetabolicgenesnoatp(konoatp==2568,:)';
nifDt=lnmetabolicgenesnoatp(konoatp==2586,:)';
nifHt=lnmetabolicgenesnoatp(konoatp==2588,:)';
nifKt=lnmetabolicgenesnoatp(konoatp==2591,:)';
petAt=lnmetabolicgenesnoatp(konoatp==2634.1,:)';
petBt=lnmetabolicgenesnoatp(konoatp==2635.1,:)';
petCt=lnmetabolicgenesnoatp(konoatp==2636.1,:)';
petDt=lnmetabolicgenesnoatp(konoatp==2637.1,:)';
petEt=lnmetabolicgenesnoatp(konoatp==2638.1,:)';
petFt=lnmetabolicgenesnoatp(konoatp==2639.1,:)';
petGt=lnmetabolicgenesnoatp(konoatp==2640.1,:)';
petHt=lnmetabolicgenesnoatp(konoatp==2641.1,:)';
psaAt=lnmetabolicgenesnoatp(konoatp==2689.1,:)';
psaBt=lnmetabolicgenesnoatp(konoatp==2690.1,:)';
psaCt=lnmetabolicgenesnoatp(konoatp==2691.1,:)';
psaDt=lnmetabolicgenesnoatp(konoatp==2692.1,:)';
psaEt=lnmetabolicgenesnoatp(konoatp==2693.1,:)';
psaFt=lnmetabolicgenesnoatp(konoatp==2694.1,:)';
psbAt=lnmetabolicgenesnoatp(konoatp==2703.1,:)';
psbBt=lnmetabolicgenesnoatp(konoatp==2704.1,:)';
psbCt=lnmetabolicgenesnoatp(konoatp==2705.1,:)';
psbDt=lnmetabolicgenesnoatp(konoatp==2706.1,:)';
psbEt=lnmetabolicgenesnoatp(konoatp==2707.1,:)';
psbFt=lnmetabolicgenesnoatp(konoatp==2708.1,:)';
psbOt=lnmetabolicgenesnoatp(konoatp==2716.1,:)';
petA2t=lnmetabolicgenesnoatp(konoatp==2634.2,:)';
petB2t=lnmetabolicgenesnoatp(konoatp==2635.2,:)';
petC2t=lnmetabolicgenesnoatp(konoatp==2636.2,:)';
petD2t=lnmetabolicgenesnoatp(konoatp==2637.2,:)';
petE2t=lnmetabolicgenesnoatp(konoatp==2638.2,:)';
petF2t=lnmetabolicgenesnoatp(konoatp==2639.2,:)';
petG2t=lnmetabolicgenesnoatp(konoatp==2640.2,:)';
petH2t=lnmetabolicgenesnoatp(konoatp==2641.2,:)';
psaA2t=lnmetabolicgenesnoatp(konoatp==2689.2,:)';
psaB2t=lnmetabolicgenesnoatp(konoatp==2690.2,:)';
psaC2t=lnmetabolicgenesnoatp(konoatp==2691.2,:)';
psaD2t=lnmetabolicgenesnoatp(konoatp==2692.2,:)';
psaE2t=lnmetabolicgenesnoatp(konoatp==2693.2,:)';
psaF2t=lnmetabolicgenesnoatp(konoatp==2694.2,:)';
psbA2t=lnmetabolicgenesnoatp(konoatp==2703.2,:)';
psbB2t=lnmetabolicgenesnoatp(konoatp==2704.2,:)';
psbC2t=lnmetabolicgenesnoatp(konoatp==2705.2,:)';
psbD2t=lnmetabolicgenesnoatp(konoatp==2706.2,:)';
psbE2t=lnmetabolicgenesnoatp(konoatp==2707.2,:)';
psbF2t=lnmetabolicgenesnoatp(konoatp==2708.2,:)';
psbO2t=lnmetabolicgenesnoatp(konoatp==2716.2,:)';
nrfAt=lnmetabolicgenesnoatp(konoatp==3385,:)';
hdrA2t=lnmetabolicgenesnoatp(konoatp==3388,:)';
hdrB2t=lnmetabolicgenesnoatp(konoatp==3389,:)';
hdrC2t=lnmetabolicgenesnoatp(konoatp==3390,:)';
norBt=lnmetabolicgenesnoatp(konoatp==4561,:)';
sorAt=lnmetabolicgenesnoatp(konoatp==5301,:)';
petJt=lnmetabolicgenesnoatp(konoatp==8906.1,:)';
pufLt=lnmetabolicgenesnoatp(konoatp==8928,:)';
pufMt=lnmetabolicgenesnoatp(konoatp==8929,:)';
NRt=lnmetabolicgenesnoatp(konoatp==10534,:)';
HAOt=lnmetabolicgenesnoatp(konoatp==10535.1,:)';
pmoAamoAt=lnmetabolicgenesnoatp(konoatp==10944,:)';
amoAt=lnmetabolicgenesnoatp(konoatp==10944.1,:)';
pmoAt=lnmetabolicgenesnoatp(konoatp==10944.2,:)';
amoBt=lnmetabolicgenesnoatp(konoatp==10945.1,:)';
pmoBt=lnmetabolicgenesnoatp(konoatp==10945.2,:)';
amoCt=lnmetabolicgenesnoatp(konoatp==10946.1,:)';
pmoCt=lnmetabolicgenesnoatp(konoatp==10946.2,:)';
dsrAoxt=lnmetabolicgenesnoatp(konoatp==11180.1,:)';
dsrAredt=lnmetabolicgenesnoatp(konoatp==11180.2,:)';
dsrBoxt=lnmetabolicgenesnoatp(konoatp==11181.1,:)';
dsrBredt=lnmetabolicgenesnoatp(konoatp==11181.2,:)';

nrfHt=lnmetabolicgenesnoatp(konoatp==15876,:)';
mxaCt=lnmetabolicgenesnoatp(konoatp==16257,:)';
mxaLt=lnmetabolicgenesnoatp(konoatp==16259,:)';
soxAt=lnmetabolicgenesnoatp(konoatp==17222,:)';
soxXt=lnmetabolicgenesnoatp(konoatp==17223,:)';
soxBt=lnmetabolicgenesnoatp(konoatp==17224,:)';
soxCt=lnmetabolicgenesnoatp(konoatp==17225,:)';
soxYt=lnmetabolicgenesnoatp(konoatp==17226,:)';
soxZt=lnmetabolicgenesnoatp(konoatp==17227,:)';
fccBt=lnmetabolicgenesnoatp(konoatp==17229,:)';
fccAt=lnmetabolicgenesnoatp(konoatp==17230,:)';
soxDt=lnmetabolicgenesnoatp(konoatp==22622,:)';
xoxFt=lnmetabolicgenesnoatp(konoatp==23995,:)';

%also fix nosZ taxa
nosZ_taxa=readmatrix("nosZ_gamma_bacteroidetes.txt");
%These are in RPKM

%add a pseudocount to 0 values
nosZ_taxa(nosZ_taxa(:)==0)=0.01;
nosZbactert=log(nosZ_taxa(1,:)./meanhousekeeping)';
nosZgammat=log(nosZ_taxa(2,:)./meanhousekeeping)';

names=['  nirB  '
       '  nirD  '
       '  nirA  '
       '  narB  '
       '  nirK  '
       '  narG  '
       '  narH  '
       '  nasC  '
       '  narI  '
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
       '  amoB  '
       '  pmoB  '
       '  amoC  '
       '  pmoC  '
       ' dsrAox '
       ' dsrAred'
       ' dsrBox '
       ' dsrBred'
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
   
%writematrix(lnmetabolicgenesnoatp, "lnmetabolicgenesnoatp.txt");
%writematrix(konoatp',"konoatp.txt");
%writematrix(uat, "uat.txt");
%writematrix(uct, "uct.txt");
%writematrix(lnhousekeepinggenes,"lnhousekeepinggenes.txt");
%writematrix(kohousekeeping,"kohousekeeping.txt");
%writematrix(lnmetabolicgenesnoatpanom, "lnmetabolicgenesnoatpanom.txt");

%Download these files to:
%/Users/sarahpreheim/Documents/Solexa_dir/Chesapeake_Bay/metagenomics/metagenomic_assembly/Mainstem/normalized_cleaned_RPKM/TPM_merged_by_KO_final_dir
%Run the following:
%perl make_dataXY_correl_matlab_umt.pl uat.txt lnmetabolicgenesnoatp.txt konoatp.txt test_uat
%upload them again and run
%dataX = table2array(test_uat_X(:,2:end));
%dataY = table2array(test_uat_Y(:,2:end));
%[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX',dataY');
%writematrix(pval, "pval_uat.txt");
%writematrix(corr_obs, "corr_obs_uat.txt");
%also do the same with the metadata to get correlations with environmental
%variables
%in:
%/Users/sarahpreheim/Documents/Solexa_dir/Chesapeake_Bay/metagenomics/metagenomic_assembly/Mainstem/normalized_cleaned_RPKM/TPM_merged_by_KO_final_dir
%Run:
%perl make_dataXY_correl_matlab_envars.pl uat.txt MATLAB_metadata_corr.txt MATLAB_names.txt test_uat_envars
%Import back into matlab
%dataX = table2array(test_uat_envars_X(:,2:end));
%dataY = table2array(test_uat_envars_Y(:,2:end));
%[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX',dataY');
%writematrix(pval, "pval_envars.txt")
%writematrix(corr_obs, "corr_obs_envars.txt")

%Read metadata for other scripts
%Some scripts use just the environmental variables
MATLAB_metadata_2=readtable("MATLAB_metadata_2.txt");
%Other scripts use all metadata data (i.e. Yearday, station etc)
metadata=readtable("MATLAB_metadata.txt", "ReadRowNames",true);
