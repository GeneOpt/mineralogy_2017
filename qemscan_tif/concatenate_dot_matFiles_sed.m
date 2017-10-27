%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder
clear all
close all
run mineral_colors.m;
    
combMin = {'Plag_A','Biot_A','Chlor_A','Ill_smec_A','Dol_A','Ilm_A'}
combMinP = {'Plag A','Biot A','Chlor A','Ill smec A','Dol A','Ilm A'}

% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\';
[nms] = dir([folder]);
matNm = {nms.name};

% for H = 3:28  % for the rock
for Hh = 5:23    % for the sediment  going to have to continue this at Hh = 5
    
    fname = dir([folder matNm{Hh} '\*.mat'])
    fileName = {fname.name}

    MAM_c = zeros(39,39,4);
    MAM_all_c = zeros(39,39);
    D_c = [];
    Ar_c = [];
    mnrlMtx_c = [];
    
    
    mn = matNm{Hh};
    mn5 = mn(1:5);
    load([folder 'concatenated files\' mn5 '.mat'])
    ptclH = struct;
    
        

%     for M = 2:length(mins)
        for M = 1:45
        islesizeM = [];
        ptclSizeM = [];
        ptclH_M = []
        
        
        
        for F = 1:length(fileName)
           
            load([folder matNm{Hh} '\' fileName{F}])
            fieldN = fieldnames(I_mtxM_all);
%         MAM_c = [MAM_c+MAM]; % here you add together the mineral association matrix (MAM)
%         MAM_all_c = [MAM_all_c+MAM_all]; % and the mineral association mtx for all sizes
%         D_c = [D_c,D]; % and concatenate the size vectors
%         Ar_c = [Ar_c,Ar]; 
%         mnrlMtx_c = [mnrlMtx_c;mnrlMtx];

            [islSize, ptclSize,ptcl_H] = islandWithin(I_mtx,I_mtxM_all.(fieldN{M}),mnrlMtx);
            islesizeM = [islesizeM,islSize];
            ptclSizeM = [ptclSizeM,ptclSize];
            ptclH_M = [ptclH_M,ptcl_H];
%             subM = MAbb.(mins_abbv{MAbb_i(M)});
%             
%             for MM = 1:length((subM))
%                 
%                 dummyMtx = I_mtxM_all.(mins{subM(MM)});
% 
%                 if MM == 1
%                     imtxLoop = dummyMtx;
%                     mDum = max(max(dummyMtx));
%                 else
%                     dummyMtx(dummyMtx~=0) = dummyMtx(dummyMtx~=0)+mDum;
%                     imtxLoop = imtxLoop+(dummyMtx);
%                     mDum = mDum + max(max(dummyMtx));
%                 end
%             end
%             I_mtxM_all_abv.(mins_abbv{MAbb_i(M)}) = imtxLoop;

        end

%         save([folder matNm{H} '\' fileName{F}],'I_mtxM_all_abv','-append')    
        isleD.(fieldN{M}) = islesizeM;
        ptclD.(fieldN{M}) = ptclSizeM;
        ptclH.(fieldN{M}) = ptclH_M;


    end

    % for the seds

    save([folder 'concatenated files\' mn5 '.mat'],'isleD','ptclD','ptclH','-append');

%     for the rock
%     save([folder 'concatenated files\' matNm{H}],'ptclD_abbv','isleD_abbv','-append');
    
end

    


