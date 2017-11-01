clear all
close all
run loadSample_specs
run mineral_colors.m
    
xValR = [16 1 2 6 7 7 8 17 3 3 4 4 9 9 18 18 12 12 13 19 19 20 20 14 14 15]
xValS = [16 1 5 2 6 7 8 17 3 4 9 10 18 11 12 13 19 20 14 15 16]   
s2rA = [1 NaN; 2 NaN; NaN NaN; 3 NaN; 4 NaN; 5 6; 7 NaN; 8 NaN; 9 10; 11 12;...
    13 14; NaN NaN; 15 16; 16 NaN; 17 18; 19 NaN; 20 21; 22 23; 24 25; 26 NaN]
% get all the folder names (each folder is a rock sample)

fOut = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\fitGSD\islands\sed'

folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
[nmsR] = dir([folderR '\*.mat']);
matNmR = {nmsR.name}

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nmsS] = dir([folderS '\*.mat']);
matNmS = {nmsS.name}

dummyV = load([folderR 'GL01RS01.mat']);
ID = dummyV.isleD;
fID = fields(ID)

colR = [0.3 0.3 0.3;1 0 0];
colS = [0.7 0.7 0.7;1 0 1];
minP = 0.7420;
maxP = 62.0;
nBins = 10;

for i = 1:length(matNmS)

    mns = matNmS{i};
    varS = load([folderS mns]);
    varB = load([folderS matNmS{i}]);
    
    mnrl_S = varS.mnrlMtx_c;
    [mnrlS,minsN,minNFull] = abbvMins(mnrl_S,5);
    islD = varS.isleD;
    
    fOutG = [fOut '\' matNmS{i}(1:5)];
    mkdir(fOutG)
   

    for M = 1:length(fID)
        
        cnt = 1;
        elS = [];

%         for k = 1:size(mnrlS,1)
%             a = mnrlS(k,:);
%             aN = a(M)./sum(a(1:38));
%             if aN>0
%                 elS(cnt) = k;
%                 cnt = cnt+1;
%             end
% 
%         end   

%         gS = varS.D_c;
        islM = (islD.(fID{M}));
        x = islM;
%         x = gS;
        x(x==0) = [];
        mA(M) = mean(x);
        stdA(M) = std(x);
        kurtA(M) = kurtosis(x);
        skewA(M) = skewness(x);

    end
 
    sampStats= [mA',stdA',kurtA',skewA'];
    return
end

% save([fOut '\paramM_GSD_all.mat'],'paramM')
    










    
    