run loadSample_specs

varS = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\GL 01.mat')
varR = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\GL01RS01.mat')
varB = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\GL01B.mat')

MAM_cS = varS.MAM_c;
MAM_cB = varB.MAM_c;
MAM_cR = varR.MAM_c;
MAM_all_cS = varS.MAM_all_c;
MAM_all_cR = varR.MAM_all_c;
MAM_all_cB = varB.MAM_all_c;

AbM = 0

if AbM == 1
    [MAM_cS,minsN,minNFull] = abbvMins(MAM_cS,4)
    [MAM_cR,minsN,minNFull] = abbvMins(MAM_cR,4)
    [MAM_cB,minsN,minNFull] = abbvMins(MAM_cB,4)
    [MAM_all_cS] = abbvMins(MAM_all_cS,3)
    [MAM_all_cR] = abbvMins(MAM_all_cR,3)
    [MAM_all_cB] = abbvMins(MAM_all_cB,3)
end





