%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder

run mineral_colors.m;
    
% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\mineral_clusters\rock\';
nmsFol = dir(folder);
matNmFl = {nmsFol.name};
matNmFl = matNmFl(3:end);

% grab the folder and look inside
for F = 1:length(matNmFl)
    inFold = [folder matNmFl{F}];
    [nms] = dir([inFold '\*.mat']);
    matNm = {nms.name};
    
    for m = 2:length(mins)

        minList = [];
        for N = 1:length(matNm)

            fname = [inFold '\' matNm{N}];
            load(fname);
            minList = [minList;clustM.(mins{m})];
            
        end
        
        clustM_c.(mins{m}) = minList;
        
    end

    mn = matNm{1};
    mn5 = mn(1:4);
    save([inFold '.mat'],...
        'clustM_c');

end

