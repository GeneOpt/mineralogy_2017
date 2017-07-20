% this function finds grains of a given mineralogy, and then computed 
% their homogeneity and returns the H index and associated grain size. If
% you want to know the H of all minerlas combined, then you ca ngo to one
% of the lower section in the script "plot_multiMineralic" from which the
% function herein is called

function [H, sz] = homog_by_mineral(MnrlMtx,elMin,limMinAb,nRow)

H = [];
sz = [];
count = 1;

for i = 1:nRow
    a = MnrlMtx(i,:);
    na = (a/sum(a));
    if na(elMin)>limMinAb      
        grain = sort(na,'descend');
        t = zeros(1,length(na)-1);
        t(1) = grain(1);
        for j = 2:length(na)-1
            t(j) = t(j-1)-(grain(j)*grain(j+1));
        end
        H(count) = t(end);
        sz(count) = sum(a);
        count = count+1;
    end
end