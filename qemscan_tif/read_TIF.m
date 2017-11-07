clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\revisedColors\sediment\MI8118-APR16\1 - 25 um\';
img_num = '';
img = 'GL 09_01';
f_img = [img img_num];
fname = [folder img img_num '.TIF'];
[I,cmap] = imread(fname);

Irgb = ind2rgb(I,cmap);
figure('position',[0.0010    0.7610    1.9200    0.9673]*1e3)
imshow(Irgb)

Ir = Irgb(:,:,1);
Ig = Irgb(:,:,2);
Ib = Irgb(:,:,3);

%% Plot select minerals

mineral = 'Ab';
minRGB = min_col(find(strcmp(mins,mineral)),:);
cols = [0 0 0];
[newIo] = keepOnly(Ir,Ig,Ib,minRGB,0,cols);
newI = newIo;
c = 0;

% mineral = 'Qtz'
% minRGB = min_col(find(strcmp(mins,mineral)),:)
% cols = []
% [newIo] = keepOnly(Ir,Ig,Ib,minRGB,0,cols);
% newI = newI+newIo;
% c = c+1
 
figure('position',[0.0010    0.7610    1.9200    0.9673]*1e3)
imshow(newI-c)



