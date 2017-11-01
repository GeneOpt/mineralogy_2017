function [n,yF,binC,m,b] = fitSL(x)
%%
minX = min(x);
maxX = max(x);
bins = 8;
binLims = linspace(minX,maxX,bins+1);
dBin = binLims(2)-binLims(1);
binC = binLims(1:end-1)+dBin/2;

yBin = [];
for i = 1:bins

    inBin = (x<binLims(i+1))&(x>binLims(i));
    yBin(i) = length(x(inBin));   

end
n = yBin/sum(yBin);
eln0 = find(n==0);

if isempty(eln0)==0
    eln0 = eln0(1);
    n = n(1:eln0-1);
    binC = binC(1:eln0-1);
end

[m,b,~,~,~] = linReg(binC,log(n));
yy = m*binC+b;
yF = exp(yy);

% figure
% plot(log2(binC/1000),(n),'k-o')
% hold on
% plot(log2(binC/1000),yF,'r-o')
% keyboard




