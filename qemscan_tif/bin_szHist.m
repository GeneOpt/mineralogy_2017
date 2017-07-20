% this sunction is called from plot_multiMineralic and takes the sizes and
% homogeneity index and bins them in a histogram fashion

function [binC,yBin] = bin_szHist(x,bins,y)
    lsa = log(x);
    mxsa = max(lsa);
    mxsa = log(976);
    mnsa = 0;
    delS = (mxsa-mnsa)/bins;
    binLims = [mnsa:delS:mxsa];
    binC = binLims(1:end-1)+(delS/2);

    for ii = 1:bins

        inBin = (lsa<binLims(ii+1))&(lsa>binLims(ii));
        yBin(ii) = length(y(inBin)); 
%         keyboard
    end
end