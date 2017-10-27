%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder
clear all
run mineral_colors.m;
scf = 25/38;
    
% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nms] = dir([folder '\*.mat']);
matNm = {nms.name}

for F = 1:length(matNm)

    load([folder matNm{F}]);
    isleD = rmfield(isleD,'Ilmc');
    ptclD = rmfield(ptclD,'Ilmc');
    ptclH = rmfield(ptclH,'Ilmc');
    save([folder matNm{F}],'isleD','ptclD','ptclH','-append')

%     mnrlM = vars.mnrlMtx_c;
%     [mnrlM_ab,minsN_abbv,minNFull_abbv] = abbvMins(mnrlM,2);
%     H = homog_by_mineral(mnrlM);
%     H_ab = homog_by_mineral(mnrlM_ab);
%     numPix = sum(mnrlM(:,2:end),2);
%     D_c = sqrt(numPix*(scf^2)*4/pi);
%     Ar_c = sqrt(numPix*(scf^2));
%     save([folder matNm{F}],'D_c','Ar_c','-append');

end

