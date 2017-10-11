% this sunction is called from plot_multiMineralic and takes the sizes and
% homogeneity index and bins them in a histogram fashion

function [binC,yBin] = bin_szMean(x,y,bins)
%     lsa = log(x);
    lsa = x
%     mxsa = max(lsa);
%     mxsa = log(976);
    mxsa = 31;
    mnsa = 0;
    delS = (mxsa-mnsa)/bins;
    binLims = [mnsa:delS:mxsa];
    binC = binLims(1:end-1)+(delS/2);

    for i = 1:bins

        inBin = (lsa<binLims(i+1))&(lsa>binLims(i));
        yBin(i) = mean(y(inBin));   

    end
end