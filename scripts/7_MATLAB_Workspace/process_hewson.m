%process_hewson_rna
load HewsonTPMfinaltable.mat
%% 
RNA(:,1)=mean(HewsonTPMfinaltable{:,41:42}')'; mon_rna(1)=10;day_rna(1)=18;year_rna(1)=2010;depth_rna(1)=13;
RNA(:,2)=mean(HewsonTPMfinaltable{:,45:48}')';mon_rna(2)=4;day_rna(2)=18;year_rna(2)=2011;depth_rna(2)=13.5;
RNA(:,3)=mean(HewsonTPMfinaltable{:,25:26}')';mon_rna(3)=5;day_rna(3)=17;year_rna(3)=2010;depth_rna(3)=13;
RNA(:,4)=mean(HewsonTPMfinaltable{:,27:28}')';mon_rna(4)=5;day_rna(4)=17;year_rna(4)=2010;depth_rna(4)=3;
RNA(:,5)=mean(HewsonTPMfinaltable{:,[43:44 49:50]}')';mon_rna(5)=5;day_rna(5)=24;year_rna(5)=2011;depth_rna(5)=18;
RNA(:,6)=mean(HewsonTPMfinaltable{:,[6:7 52] }')';mon_rna(6)=6;day_rna(6)=14;year_rna(6)=2011;depth_rna(6)=3;
RNA(:,7)=mean(HewsonTPMfinaltable{:,53:55}')';mon_rna(7)=6;day_rna(7)=14;year_rna(7)=2011;depth_rna(7)=18;
RNA(:,8)=mean(HewsonTPMfinaltable{:,29:30}')';mon_rna(8)=6;day_rna(8)=7;year_rna(8)=2010;depth_rna(8)=17.5;
RNA(:,9)=mean(HewsonTPMfinaltable{:,8:10}')';mon_rna(9)=7;day_rna(9)=8;year_rna(9)=2011;depth_rna(9)=3;
RNA(:,10)=mean(HewsonTPMfinaltable{:,31:32}')';mon_rna(10)=7;day_rna(10)=11;year_rna(10)=2010;depth_rna(10)=17.5;
RNA(:,11)=mean(HewsonTPMfinaltable{:,33:34}')';mon_rna(11)=7;day_rna(11)=11;year_rna(11)=2010;depth_rna(11)=3;
RNA(:,12)=mean(HewsonTPMfinaltable{:,11:12}')';mon_rna(12)=8;day_rna(12)=30;year_rna(12)=2011;depth_rna(12)=18;
RNA(:,13)=mean(HewsonTPMfinaltable{:,13:16}')';mon_rna(13)=8;day_rna(13)=30;year_rna(13)=2011;depth_rna(13)=3;
RNA(:,14)=mean(HewsonTPMfinaltable{:,17:20}')';mon_rna(14)=8;day_rna(14)=8;year_rna(14)=2011;depth_rna(14)=18;
RNA(:,15)=mean(HewsonTPMfinaltable{:,39:40}')';mon_rna(15)=8;day_rna(15)=30;year_rna(15)=2010;depth_rna(15)=20;
RNA(:,16)=mean(HewsonTPMfinaltable{:,37:38}')';mon_rna(16)=8;day_rna(16)=30;year_rna(16)=2010;depth_rna(16)=3;
RNA(:,17)=mean(HewsonTPMfinaltable{:,35:36}')';mon_rna(17)=8;day_rna(17)=5;year_rna(17)=2010;depth_rna(17)=22;
RNA(:,18)=mean(HewsonTPMfinaltable{:,21:24}')';mon_rna(18)=9;day_rna(18)=21;year_rna(18)=2011;depth_rna(18)=17;
Yearday_rna=datenum(year_rna',mon_rna',day_rna')-datenum(year_rna',1,1)+1;
%% 
%
hydro_samp_rna=[90	91 8 2 105 108 114 16 124 40 34 200 194 192 82 76 74 208]';
    

%% 
nosZrna=RNA(HewsonTPMfinaltable{:,1}==376,:);
nifDrna=RNA(HewsonTPMfinaltable{:,1}==2586,:);
nifHrna=RNA(HewsonTPMfinaltable{:,1}==2588,:);
nifKrna=RNA(HewsonTPMfinaltable{:,1}==2591,:);
psaBrna=RNA(HewsonTPMfinaltable{:,1}==2690,:);
haorna=RNA(HewsonTPMfinaltable{:,1}==10535,:);
amoArna=RNA(HewsonTPMfinaltable{:,1}==10944,:);
amoBrna=RNA(HewsonTPMfinaltable{:,1}==10945,:);
dsrArna=RNA(HewsonTPMfinaltable{:,1}==11180,:);
dsrBrna=RNA(HewsonTPMfinaltable{:,1}==11180,:);
%% 
amoABrna=amoArna+amoBrna;
dsrABrna=dsrArna+dsrBrna;
nifDHKrna=nifDrna+nifHrna+nifKrna;
figure(12)
clf;
subplot(231)
scatter(Yearday_rna,-depth_rna,40,log10(psaBrna))
hold on
scatter(Yearday_rna(psaBrna>0),-depth_rna(psaBrna>0),40,log10(psaBrna(psaBrna>0)),'filled');colorbar
title('(a) log_{10} {\it psaB} RNA, CB4.3')
xlabel('Time in decimal day')
ylabel('Depth in m')
subplot(232)
scatter(Yearday_rna,-depth_rna,40,log10(haorna))
hold on
scatter(Yearday_rna(haorna>0),-depth_rna(haorna>0),40,log10(haorna(haorna>0)),'filled');colorbar
title('(b) log_{10} {\it hao} RNA, CB4.3')
xlabel('Time in decimal day')
ylabel('Depth in m')
subplot(233)
scatter(Yearday_rna,-depth_rna,40,log10(amoABrna))
hold on
scatter(Yearday_rna(amoABrna>0),-depth_rna(amoABrna>0),40,log10(amoABrna(amoABrna>0)),'filled');colorbar
title('(c) log_{10} {\it amoAB} RNA, CB4.3')
xlabel('Time in decimal day')
ylabel('Depth in m')
subplot(234)
scatter(Yearday_rna,-depth_rna,40,log10(nosZrna))
hold on
scatter(Yearday_rna(nosZrna>0),-depth_rna(nosZrna>0),40,log10(nosZrna(nosZrna>0)),"filled");colorbar
title('(d) log_{10} {\it nosZ} RNA, CB4.3')
xlabel('Time in decimal day')
ylabel('Depth in m')
subplot(235)
scatter(Yearday_rna,-depth_rna,40,log10(dsrABrna))
hold on
scatter(Yearday_rna(dsrABrna>0),-depth_rna(dsrABrna>0),40,log10(dsrABrna(dsrABrna>0)),"filled");colorbar
title('(e) log_{10} {\it dsrAB} RNA, CB4.3')
xlabel('Time in decimal day')
ylabel('Depth in m')
subplot(236)
scatter(Yearday_rna,-depth_rna,40,log10(nifDHKrna))
hold on
scatter(Yearday_rna(nifDHKrna>0),-depth_rna(nifDHKrna>0),40,log10(nifDHKrna(nifDHKrna>0)),"filled");colorbar
title('(f) log_{10} {\it nifDHK} RNA, CB4.3')
xlabel('Time in decimal day')
ylabel('Depth in m')




