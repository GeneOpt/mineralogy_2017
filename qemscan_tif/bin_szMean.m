% this sunction is called from plot_multiMineralic and takes the sizes and
% homogeneity index and bins them in a histogram fashion

function [binC,yBin] = bin_szMean(sz_Arr,T_Arr,bins)
    lsa = log(sz_Arr);
    mxsa = max(lsa);
    mxsa = log(976);
    mnsa = 0;
    delS = (mxsa-mnsa)/bins;
    binLims = [mnsa:delS:mxsa];
    binC = binLims(1:end-1)+(delS/2);

    for i = 1:bins

        inBin = (lsa<binLims(i+1))&(lsa>binLims(i));
        yBin(i) = mean(T_Arr(inBin));   

    end
end