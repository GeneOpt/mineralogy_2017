function [] = testRGB(I_mtx,count)

R = ones(size(I_mtx));
G = ones(size(I_mtx));
B = ones(size(I_mtx));

R(I_mtx==count)=0;
G(I_mtx==count)=0;
B(I_mtx==count)=1;

RGB = zeros(size(I_mtx,1),size(I_mtx,2),3);
RGB(:,:,1) = R;
RGB(:,:,2) = G;
RGB(:,:,3) = B;

figure
imshow(RGB)