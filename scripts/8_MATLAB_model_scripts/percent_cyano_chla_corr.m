

cyanophoto(1,:)=nlpetAt1;
cyanophoto(2,:)=nlpetBt1;
cyanophoto(3,:)=nlpetCt1;
cyanophoto(4,:)=nlpetDt1;
cyanophoto(5,:)=nlpetEt1;
cyanophoto(6,:)=nlpetFt1;
cyanophoto(7,:)=nlpetGt1;
cyanophoto(8,:)=nlpetHt1;
cyanophoto(9,:)=nlpetJt1;

cyanophoto(10,:)=nlpsaAt1;
cyanophoto(11,:)=nlpsaBt1;
cyanophoto(12,:)=nlpsaCt1;
cyanophoto(13,:)=nlpsaDt1;
cyanophoto(14,:)=nlpsaEt1;
cyanophoto(15,:)=nlpsaFt1;

cyanophoto(16,:)=nlpsbAt1;
cyanophoto(17,:)=nlpsbBt1;
cyanophoto(18,:)=nlpsbCt1;
cyanophoto(19,:)=nlpsbDt1;
cyanophoto(20,:)=nlpsbEt1;
cyanophoto(21,:)=nlpsbFt1;
cyanophoto(22,:)=nlpsbOt1;

eukphoto(1,:)=nlpetAt2;
eukphoto(2,:)=nlpetBt2;
eukphoto(3,:)=nlpetCt2;
eukphoto(4,:)=nlpetDt2;
eukphoto(5,:)=nlpetEt2;
eukphoto(6,:)=nlpetFt2;
eukphoto(7,:)=nlpetGt2;
eukphoto(8,:)=nlpetHt2;
eukphoto(9,:)=nlpetJt2;

eukphoto(10,:)=nlpsaAt2;
eukphoto(11,:)=nlpsaBt2;
eukphoto(12,:)=nlpsaCt2;
eukphoto(13,:)=nlpsaDt2;
eukphoto(14,:)=nlpsaEt2;
eukphoto(15,:)=nlpsaFt2;

eukphoto(16,:)=nlpsbAt2;
eukphoto(17,:)=nlpsbBt2;
eukphoto(18,:)=nlpsbCt2;
eukphoto(19,:)=nlpsbDt2;
eukphoto(20,:)=nlpsbEt2;
eukphoto(21,:)=nlpsbFt2;
eukphoto(22,:)=nlpsbOt2;



percentphoto=mean(cyanophoto./(cyanophoto+eukphoto));

answer=ttest(percentphoto(pval_chla(1:22)>0.05),percentphoto(pval_chla(1:22)<0.05));
%This is not significant, so it doesn't seem like the large taxa breakdown
%explains this