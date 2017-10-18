%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder
clear all
run mineral_colors.m;
    
% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
[nms] = dir([folder '\*.mat']);
matNm = {nms.name}

for F = 1:length(matNm)
    
    vars = load([folder matNm{F}]);
    mnrlM = vars.mnrlMtx_c;
    [mnrlM_ab,minsN_abbv,minNFull_abbv] = abbvMins(mnrlM,2);
    H = homog_by_mineral(mnrlM);
    H_ab = homog_by_mineral(mnrlM_ab);
    save('[folder matNm{F}]','H','H_ab','-append')
    
end

