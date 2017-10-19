%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder
clear all
close all
run mineral_colors.m;
    
% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\';
[nms] = dir([folder]);
matNm = {nms.name};

% for H = 3:28  % for the rock
for H = 3:23    % for the sediment
    
    fname = dir([folder matNm{H} '\*.mat'])
    fileName = {fname.name}

    MAM_c = zeros(39,39,4);
    MAM_all_c = zeros(39,39);
    D_c = [];
    Ar_c = [];
    mnrlMtx_c = [];
    

    for M = 2:length(mins)
        islesizeM = [];
        ptclSizeM = [];
        for F = 1:length(fileName)
        
            load([folder matNm{H} '\' fileName{F}])
        
%         MAM_c = [MAM_c+MAM]; % here you add together the mineral association matrix (MAM)
%         MAM_all_c = [MAM_all_c+MAM_all]; % and the mineral association mtx for all sizes
%         D_c = [D_c,D]; % and concatenate the size vectors
%         Ar_c = [Ar_c,Ar]; 
%         mnrlMtx_c = [mnrlMtx_c;mnrlMtx];
  
            [islSize, ptclSize] = islandWithin(I_mtx,I_mtxM_all.(mins{M}));
            islesizeM = [islesizeM,islSize];
            ptclSizeM = [ptclSizeM,ptclSize];
            
        end
            
        isleD.(mins{M}) = islesizeM;
        ptclD.(mins{M}) = ptclSizeM;
        figure
        plot(ptclSizeM,islesizeM,'.')
        grid on
        xlabel('ptclSizeM')
        ylabel('isleSizeM')
        title(mins{M})


    end

    % for the seds
    mn = matNm{H};
    mn5 = mn(1:5);
    save([folder 'concatenated files\' mn5 '.mat'],'','-append');

    % for the rock
%     save([folder 'concatenated files\' matNm{H}],'MAM_c','MAM_all_c','D_c','Ar_c','mnrlMtx_c');
    
end

    
%     for m = 2:length(mins)
% 
%         minList = [];
%         for N = 1:length(matNm)
% 
%             fname = [inFold '\' matNm{N}];
%             load(fname);
%             minList = [minList;clustM.(mins{m})];
%             
%         end
%         
%         clustM_c.(mins{m}) = minList;
%         
%     end


