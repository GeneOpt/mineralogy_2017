function[rs] = rsquare(yi, f)
yibar = sum(yi) / length(yi);

SStot = sum((yi - yibar).^2);
SSres = sum((yi - f).^2);

rs = 1 - (SSres/SStot);