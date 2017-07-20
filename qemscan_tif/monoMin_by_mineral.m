% this function finds grains that are more than a given fraction (like
% 0.95) monomineralic. To normalize, the histogram needs to be run on all
% minerals of any composition, then a hist of sz needs to be run, and the
% outcome of the two divided. 
function [sz] = monoMin_by_mineral(MnrlMtx,elMin,limMinAb,nRow)

sz = [];
count = 1;

for i = 1:nRow
    a = MnrlMtx(i,:);
    na = (a/sum(a));
    if na(elMin)>limMinAb      
        sz(count) = sum(a);
        count = count+1;
    end
end