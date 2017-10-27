clear all
close all
run loadSample_specs
run mineral_colors.m

%% sample specs for rock

% samplesR = {'GL01_1','GL02_2','GL04_1','GL05_1','GL06_2','GL06_4','GL07_1',...
%     'GL08_2','GL09_A','GL09_B','GL10_1','GL10_2','GL11_1','GL11_2','GL14_2',...
%     'GL14_3','GL16_1','GL16_2','GL17_1','GL18_1','GL18_2','GL19_1','GL19_2',...
%     'GL20_1','GL20_2','GL21_1'}; the entire suite
samplesR = {'GL01_1','GL02_2','GL04_1','GL05_1','GL06_2','GL06_4','GL07_1',...
    'GL08_2','GL09_A','GL09_B','GL10_1','GL10_2','GL11_1','GL11_2','GL14_2',...
    'GL14_3','GL16_1','GL16_2','GL17_1','GL18_1','GL18_2','GL19_1','GL19_2',...
    'GL20_1','GL20_2','GL21_1'};
labelR = {'1','2','4','5','6a','6b','7','8','9a','9b','10a','10b','11a',...
    '11b','14a','14b','16a','16b','17','18a','18b','19a','19b','20a','20b','21'}
xValR = [16 1 2 6 7 7 8 17 3 3 4 4 9 9 18 18 12 12 13 19 19 20 20 14 14 15]

for i = 1:length(samplesR)
    keepLabR(i) = find(strcmp(varsR(:,8),samplesR(i)))
end

stR = strcmp(varsR(keepLabR,2),'S')
nstR = strcmp(varsR(keepLabR,2),'NS')
gtR = varsR(keepLabR,2)
prR = strcmp(varsR(keepLabR,3),'P')
msR = strcmp(varsR(keepLabR,3),'MS')
groupR = datR(keepLabR-1,7)


%% sample specs for sediment
% samplesS = {'GL 01','GL 02','GL 03','GL 04','GL 05','GL 06','GL 07','GL 08',...
%     'GL 09','GL 10','GL 11','GL 13','GL 14','GL 15','GL 16','GL 17','GL 18',...
%     'GL 19','GL 20','GL 21','GL1_1'}; the entire suite
samplesS = {'GL 01','GL 02','GL 03','GL 04','GL 05','GL 06','GL 07','GL 08',...
    'GL 09','GL 10','GL 11','GL 13','GL 14','GL 15','GL 16','GL 17','GL 18',...
    'GL 19','GL 20','GL 21','GL1_1'};
labelS = {'1','2','3','4','5','6','7','8','9','10','11','13',...
    '14','15','16','17','18','19','20','21','BH_1'}
xValS = [16 1 5 2 6 7 8 17 3 4 9 10 18 11 12 13 19 20 14 15 16]   


for i = 1:length(samplesS)
    keepLabS(i) = find(strcmp(varsS(:,6),samplesS(i)))
end

stS = strcmp(varsS(keepLabS,2),'S')
nstS = strcmp(varsS(keepLabS,2),'NS')
gtS = varsS(keepLabS,2)
mxS = strcmp(varsS(keepLabS,3),'MX')
msS = strcmp(varsS(keepLabS,3),'MS')
groupS = datS(keepLabS-1,5)
labelS = labelS(keepLabS-1)
%%

folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
[nmsR] = dir([folderR '\*.mat']);
matNmR = {nmsR.name}

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nmsS] = dir([folderS '\*.mat']);
matNmS = {nmsS.name}


min1 = find(strcmp(minsN,'Chl Fe'))
min2 = find(strcmp(minsN,'Ab'))

%%
f1 = figure
for Hs = 1:length(xValS)    % for the sediment

    mns = matNmS{Hs};
    varS = load([folderS mns]);
    
    MAM_cS = varS.MAM_c;
    MAM_all_cS = varS.MAM_all_c;
    
    tArrS = MAM_all_cS(min1,:);
    numPtS = sum(tArrS);
    yVal = (MAM_all_cS(min1,min2)/numPtS*100);

    colA = [0.7 0.7 0.7;1 0.7 0.7];
    mt = {'^','o'};

    hold on
    plot(xValS(Hs),yVal,['k' mt{mxS(Hs)+1}],'markerfacecolor', colA(stS(Hs)+1,:));
    text(xValS(Hs),yVal,labelS{Hs});

end
    
for Hr = 1:length(xValR)    % for the sediment

    mnr = matNmR{Hr};
    varR = load([folderR mnr]);
    
    MAM_cR = varR.MAM_c;
    MAM_all_cR = varR.MAM_all_c;

    tArrR = MAM_all_cR(min1,:);
    numPtS = sum(tArrR);
    yVal = (MAM_all_cR(min1,min2)/numPtS*100);

    colA = [0.3 0.3 0.3;1 0.2 0.2];

    hold on
    plot(xValR(Hr),yVal,['k' mt{msR(Hr)+1}],'markerfacecolor', colA(stR(Hr)+1,:))
    text(xValR(Hr),yVal,labelR{Hr})

end
    
grid on
title(['Percent ' minsN{min1} ' neighbouring ' minsN{min2}])
ylabel('%')
set(gca,'XtickLabel',[],'fontsize',18)
saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\binary\' mins{min1} '_' mins{min2}])
savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\binary\' mins{min1} '_' mins{min2}])

