clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\rock_frags\second_set_M18182\1-25um\'
img_num = ''
img = 'GL1_12_25'
f_img = [img img_num]
fname = [folder img img_num '.TIF']
[I,cmap] = imread(fname);

Irgb = ind2rgb(I,cmap);
figure
imshow(Irgb)

Ir = Irgb(:,:,1);
Ig = Irgb(:,:,2);
Ib = Irgb(:,:,3);

% Im = im2double(I);
% [r c p] = size(Im);
% 
% imR = squeeze(Im(:,:,1));
% imG = squeeze(Im(:,:,2));
% imB = squeeze(Im(:,:,3));

imBinaryR = im2bw(Ir,graythresh(Ir));
imBinaryG = im2bw(Ir,graythresh(Ig));
imBinaryB = im2bw(Ir,graythresh(Ib));
imBinary = imcomplement(imBinaryR&imBinaryG&imBinaryB);
figure
imshow(imBinary)