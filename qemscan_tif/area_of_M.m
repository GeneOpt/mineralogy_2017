% this function finds the amount of mineral M in any grain
function [aM,sz] = area_of_M(MnrlMtx,elMin,nRow)

sz = [];

for i = 1:nRow
    a = MnrlMtx(i,:);
    aM(i) = a(elMin);
    sz(i) = sum(a);
end