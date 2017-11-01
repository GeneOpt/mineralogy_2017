clc
close all
clear variables

run mineral_colors
[datR,varsR] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',1,'B1:I27');
[datS,varsS] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',2,'B1:I27');
RSArr = [1,2,4,5,6,7,8,9,9,10,10,11,11,13,13,15,15,16,17,17,18,18,19,19,20];
SPIS =  [4,1,2,5,6,10,14,8,9,13,18,22,12,3,7,11,16,20,15,19,4];
SPIR = [4,1,5,6,10,10,14,8,9,9,13,13,18,18,12,12,7,7,11,16,16,20,20,15,15,19];


xValR = [16 1 2 6 7 7 8 17 3 3 4 4 9 9 18 18 12 12 13 19 19 20 20 14 14 15]
xValS = [16 1 5 2 6 7 8 17 3 4 9 10 18 11 12 13 19 20 14 15 16]   
s2rA = [1 NaN; 2 NaN; NaN NaN; 3 NaN; 4 NaN; 5 6; 7 NaN; 8 NaN; 9 10; 11 12;...
    13 14; NaN NaN; 15 16; 16 NaN; 17 18; 19 NaN; 20 21; 22 23; 24 25; 26 NaN]
%% sample specs for rock

% samplesR = {'GL01_1','GL02_2','GL04_1','GL05_1','GL06_2','GL06_4','GL07_1',...
%     'GL08_2','GL09_A','GL09_B','GL10_1','GL10_2','GL11_1','GL11_2','GL14_2',...
%     'GL14_3','GL16_1','GL16_2','GL17_1','GL18_1','GL18_2','GL19_1','GL19_2',...
%     'GL20_1','GL20_2','GL21_1'}; the entire suite
samplesR = {'GL01_1','GL02_2','GL04_1','GL05_1','GL06_2','GL06_4','GL07_1',...
    'GL08_2','GL09_A','GL09_B','GL10_1','GL10_2','GL11_1','GL11_2','GL14_2',...
    'GL14_3','GL16_1','GL16_2','GL17_1','GL18_1','GL18_2','GL19_1','GL19_2',...
    'GL20_1','GL20_2','GL21_1'};
labelR = {'1','2','4','5','6a','6b','7','8','9a','9b','10a','10b','11a',...
    '11b','14a','14b','16a','16b','17','18a','18b','19a','19b','20a','20b','21'}

for i = 1:length(samplesR)
    keepLabR(i) = find(strcmp(varsR(:,8),samplesR(i)))
end

stR = strcmp(varsR(keepLabR,2),'S')
nstR = strcmp(varsR(keepLabR,2),'NS')
gtR = varsR(keepLabR,2)
prR = strcmp(varsR(keepLabR,3),'P')
msR = strcmp(varsR(keepLabR,3),'MS')
groupR = datR(keepLabR-1,7)



%% sample specs for sediment
% samplesS = {'GL 01','GL 02','GL 03','GL 04','GL 05','GL 06','GL 07','GL 08',...
%     'GL 09','GL 10','GL 11','GL 13','GL 14','GL 15','GL 16','GL 17','GL 18',...
%     'GL 19','GL 20','GL 21','GL1_1'}; the entire suite
samplesS = {'GL 01','GL 02','GL 03','GL 04','GL 05','GL 06','GL 07','GL 08',...
    'GL 09','GL 10','GL 11','GL 13','GL 14','GL 15','GL 16','GL 17','GL 18',...
    'GL 19','GL 20','GL 21','GL1_1'};
labelS = {'1','2','3','4','5','6','7','8','9','10','11','13',...
    '14','15','16','17','18','19','20','21','BH_1'}

for i = 1:length(samplesS)
    keepLabS(i) = find(strcmp(varsS(:,6),samplesS(i)))
end

stS = strcmp(varsS(keepLabS,2),'S')
nstS = strcmp(varsS(keepLabS,2),'NS')
gtS = varsS(keepLabS,2)
mxS = strcmp(varsS(keepLabS,3),'MX')
msS = strcmp(varsS(keepLabS,3),'MS')
groupS = datS(keepLabS-1,5)
labelS = labelS(keepLabS-1)
