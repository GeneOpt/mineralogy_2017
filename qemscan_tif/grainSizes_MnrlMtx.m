function [szA] = grainSizes_MnrlMtx(MnrlMtx,nRow)

for i=1:nRow
    a = MnrlMtx(i,:);
    szA(i) = sum(a);
end


     