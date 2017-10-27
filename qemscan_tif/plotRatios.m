clear variables
close all

run mineral_colors

min1 = find(strcmp(minsN,'Chl Fe'))
min2 = find(strcmp(minsN,'Bti low'))

vars1 = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\GL01RS01.mat');
vars2 = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\GL 01.mat')

mm1 = vars1.mnrlMtx_c;
mm2 = vars2.mnrlMtx_c;

count = 1
for i = 1:length(mm1)
    a11 = mm1(i,min1);
    a12 = mm1(i,min2);
    aT = (a11+a12)./sum(mm1(i,:));
    if aT>0.1
        elR(count) = i;
        r1(count) = a11/a12;
        count = count+1;
    end
end

count = 1
for i = 1:length(mm2)
    a11 = mm2(i,min1);
    a12 = mm2(i,min2);
    aT = (a11+a12)./sum(mm2(i,:));
    if aT>0.1
        elS(count) = i;
        r2(count) = a11/a12;
        count = count+1;
    end
end
    


el1 = isnan(r1);
r1(el1)=[];
r1(isinf(r1)) = 1e4;
r1(r1==0) = 1e-4;
el2 = isnan(r2)
r2(el2)=[];
r2(isinf(r2)) = 1e4
r2(r2==0) = 1e-4;

D_r = vars1.D_c(elR)
D_r(el1) = [];
D_s = vars2.D_c(elS)
D_s(el2) = [];

%%
close all
figure
subplot(2,1,1)
h = histogram(log(r1),'Normalization','pdf')
title('Rock')

subplot(2,1,2)
histogram(log(r2),h.BinEdges,'Normalization','pdf')
title('Sed')

%%
figure
plot(D_r,log(r1),'r.')
hold on
plot(D_s,log(r2),'k.')
grid on


n1 = length(r1)
n2 = length(r2)
s_p = sqrt((((n1-1)*std(r1)^2)+((n2-1)*std(r2)^2))/(n1+n2-2))
t = (mean(r1)-mean(r2))/(s_p*sqrt((1/n1)+(1/n2)))

%%
[binC,binC_S,yBinR,mBinR,stdBinR,sDat] = bin_szHist(D_r,log(r1));
[binC,binC_S,yBinS,mBinS,stdBinS,sDat] = bin_szHist(D_s,log(r2));

figure
plot(binC,mBinR,'r-o')
hold on
plot(binC,mBinS,'k-^')


