%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder
clear all
run mineral_colors.m;
scf = 25/38;
    
MAbb_i = [4 6 8 11 13 18]
combMin = {'Plag_A','Biot_A','Chlor_A','Ill_smec_A','Dol_A','Ilm_A'}
combMinP = {'Plag A','Biot A','Chlor A','Ill smec A','Dol A','Ilm A'}
MAbb.(combMin{1}) = [4 5 6 7]
MAbb.(combMin{2}) = [9 10 11 12]
MAbb.(combMin{3}) = [14 15]
MAbb.(combMin{4}) = [18 19]
MAbb.(combMin{5}) = [21 22]
MAbb.(combMin{6}) = [27 28]

% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\';
[nms] = dir([folder]);
matNm = {nms.name}

% for F = 3:28 %rock
for F = 3:23 %sed
    
    fName = dir([folder matNm{F} '\*.mat'])
    fN = {fName.name}
    
    for N = 1:length(fN)
        
        fullName = [folder matNm{F} '\' fN{N}]
        vars = load(fullName);
        mtx = vars.I_mtx_C;
        for i = 1:length(MAbb_i)
            M_i = MAbb_i(i);
            fMin = MAbb.(combMin{i});
            nS = length(fMin);
            dummyMtx  = zeros(size(mtx));
            for j = 1:nS
                dummyMtx(mtx==fMin(j)) = 39+i;
            end
            I_mtx_C_abbv.(combMin{i}) = dummyMtx;
        end
        save(fullName,'I_mtx_C_abbv','-append')

    end

end

