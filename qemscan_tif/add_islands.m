%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder
clear all
close all
run mineral_colors.m;
    
% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\';
[nms] = dir([folder]);
matNm = {nms.name};

% for H = 3:28  % for the rock
for H = 3   % for the sediment
    
    fname = dir([folder matNm{H} '\*.mat'])
    fileName = {fname.name}

    MAM_c = zeros(39,39,4);
    MAM_all_c = zeros(39,39);
    D_c = [];
    Ar_c = [];
    mnrlMtx_c = [];
    
    for F = 1:length(fileName)
    
    MAbb = struct
    MAbb.(mins_abbv{4}) = [4 5 6 7]
    MAbb.(mins_abbv{6}) = [9 10 11 12]
    MAbb.(mins_abbv{8}) = [14 15]
    MAbb.(mins_abbv{11}) = [18 19]
    MAbb.(mins_abbv{13}) = [21 22]
    MAbb.(mins_abbv{18}) = [27 28]
    MAbb_i = [4 6 8 11 13 18]

        for M = 1:length(MAbb_i)
   
            load([folder matNm{H} '\' fileName{F}])

            subM = MAbb.(mins_abbv{MAbb_i(M)});
            
            for MM = 1:length((subM))
                
                dummyMtx = I_mtxM_all.(mins{subM(MM)});

                if MM == 1
                    imtxLoop = dummyMtx;
                    mDum = max(max(dummyMtx));
                else
                    dummyMtx(dummyMtx~=0) = dummyMtx(dummyMtx~=0)+mDum;
                    imtxLoop = imtxLoop+(dummyMtx);
                    mDum = mDum + max(max(dummyMtx));
                end
            end
            
            I_mtxM_all_abv.(mins_abbv{MAbb_i(M)}) = imtxLoop;

        end
return
        save([folder matNm{H} '\' fileName{F}],'I_mtxM_all_abv','-append')    

    end

    
end

    %%
clear variables
clc
run mineral_colors
load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\GL01RS01\GL01RS01.mat')
dd = I_mtxM_all_abv.Plagc;
[5 10 16]
% dd(dd~=0)=1;
% imagesc(dd)
I_mtx_C(find(I_mtx==16))
mins(ans)

return
[D_isl, D_Ptcl] = islandWithin(I_mtx,dd)
figure
plot(D_Ptcl,D_isl,'o')
refline(1,0)












