% this sunction is called from plot_multiMineralic and takes the sizes and
% homogeneity index and bins them in a histogram fashion

function [binC,binC_S,yBin,mBin,stdBin,sDat,sBin] = bin_szHist(x,y)
    nPhi = log2(x/1000);
    binLims = -[10.5:-1:5.5];
    bins = length(binLims)-1;
    delS = (binLims(2)-binLims(1));
    binC = binLims(1:end-1)+(delS/2);
    binC_S = {'fine_clay','med_clay','silt_clay','vf_silt','f_silt'};
    sDat = struct();
    
    for ii = 1:bins

        inBin = (nPhi<binLims(ii+1))& (nPhi>binLims(ii));
        yB = y(inBin);
        yBin(ii) = length(yB);
        sBin(ii) = sum(yB);
        mBin(ii) = mean(yB);
        stdBin(ii) = std(yB);
        sDat.bin{ii} = inBin;

        
    end
end