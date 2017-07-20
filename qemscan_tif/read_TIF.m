clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\rock_frags\second_set_M18182\1-25um\';
img_num = '';
img = 'GL1_12_25';
f_img = [img img_num];
fname = [folder img img_num '.TIF'];
[I,cmap] = imread(fname);

Irgb = ind2rgb(I,cmap);
figure
imshow(Irgb)

Ir = Irgb(:,:,1);
Ig = Irgb(:,:,2);
Ib = Irgb(:,:,3);

return
%% Plot select minerals
figure
mineral = 'Qtz';
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
 
figure
imshow(newI-c)



