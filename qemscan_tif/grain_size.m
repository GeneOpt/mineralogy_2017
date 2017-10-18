clc
close all
clear variables

run mineral_colors

sf = 25/38;

% here you are taking the .mat files generated from 'find_isolatedGrains',
% where minerals are already tagged and numbered
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\';
[nms] = dir([folder]);
matNm = {nms.name};

for H = 3:23
    
    fname = dir([folder matNm{H} '\*.mat'])
    fileName = {fname.name}
    
    for F = 1:length(fileName)

        load([folder matNm{H} '\' fileName{F}])
        numG = max(max(I_mtx)); % there are numG grains

            for N = 1:numG

                numPix = sum(sum(I_mtx==N));
                D(N) = sqrt(numPix*(sf^2)*4/pi);
                Ar(N) = sqrt(numPix*(sf^2));

            end
        save([folder matNm{H} '\' fileName{F}],'D','Ar','-append')
    end
end
