clc
close all
clear variables

run mineral_colors

% here you are taking the .mat files generated from 'find_isolatedGrains',
% where minerals are already tagged and numbered
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\';
[nms] = dir([folder]);
matNm = {nms.name};

for H = 3:23
    
    fname = dir([folder matNm{H} '\*.mat'])
    fileName = {fname.name}
    
    for F = 1:length(fileName)

        I_mtx = []
        I_mtx_C = []
        load([folder matNm{H} '\' fileName{F}])
        numG = max(max(I_mtx)); % there are numG grains


        %%

        minAsMtx = zeros(length(mins),length(mins),4); % this matrix will be updated by everytime two grains are touching
        SArr = [0,2.^([-9.5,-8.5,-7.5,0])*1000];
        for G = 1:numG
            if D(G)<SArr(2) && D(G)>SArr(1)
                S = 1;
            elseif D(G)<SArr(3) && D(G)>SArr(2)
                S = 2;
            elseif D(G)<SArr(4) && D(G)>SArr(3)
                S = 3;
            elseif D(G)<SArr(5) && D(G)>SArr(4)
                S = 4;
            end
            elG = find(I_mtx==G);
            rw = [];
            cl = [];
            for k = 1:length(elG) % this loop simply takes the grain tag and turns it into row col
                Mk = elG(k);
                rc = ind2rc(size(I_mtx,1),Mk);
                rw(k) = int32(rc(1));
                cl(k) = int32(rc(2)); 
            end
            minR = min(rw)-1; % these next five lines takes the mineral in the I_mtx and the zeros around it
            maxR = max(rw)+1;
            minC = min(cl)-1;
            maxC = max(cl)+1;
            IGC = I_mtx_C([minR:maxR],[minC:maxC]);
            IGC(IGC==0)=1; % here you set the zeros to be 1
            size_j = size(IGC,2); % know the size of the enclosing matrix
            size_i = size(IGC,1);
            UP_G = IGC([1:end-1],:); % here you're just cutting parts of the matrices for future shifts
            DN_G = IGC([2:end],:);
            Ri_G = IGC(:,2:end);
            Le_G = IGC(:,1:end-1);
            DLU = IGC(1:end-1,1:end-1);
            DLD = IGC(2:end,2:end);
            DRD = IGC(1:end-1,2:end);
            DRU = IGC(2:end,1:end-1);

            for ii= 1:size_i-1   % in these lines you shift the matrix up, down lefft right, or diagonal
                for jj = 1:size_j   % to get the min number from each cut matrix and updated the association mtx
                    minAsMtx(UP_G(ii,jj),DN_G(ii,jj),S) = minAsMtx(UP_G(ii,jj),DN_G(ii,jj),S)+1;
                end 
            end

            for ii= 1:size_i
                for jj = 1:size_j-1
                    minAsMtx(Ri_G(ii,jj),Le_G(ii,jj),S) = minAsMtx(Ri_G(ii,jj),Le_G(ii,jj),S)+1; 
                end 
            end

            for ii= 1:size_i-1
                for jj = 1:size_j-1
                    minAsMtx(DLU(ii,jj),DLD(ii,jj),S) = minAsMtx(DLU(ii,jj),DLD(ii,jj),S)+1; 
                    minAsMtx(DRU(ii,jj),DRD(ii,jj),S) = minAsMtx(DRU(ii,jj),DRD(ii,jj),S)+1; 
                end 
            end

        end
        %% here you have to add up the upper and lower matrices, then make it symetric
        MAM = zeros(size(minAsMtx));
        for S = 1:4
            maM = minAsMtx(:,:,S);
            minAsMtx_S = triu(maM ,1)+triu(maM',1);
            minAsMtx_S2 = minAsMtx_S+minAsMtx_S';
            MAM(:,:,S) = minAsMtx_S2+diag(diag(maM));
        end
        MAM_all = sum(MAM,3);
        save([folder matNm{H} '\' fileName{F}],'MAM','MAM_all','SArr','-append')

    end
end




