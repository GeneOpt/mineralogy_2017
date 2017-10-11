function[rs] = rsquare(yi, fi)
yibar = sum(yi) / length(yi);

SStot = sum((yi - yibar).^2);
SSres = sum((yi - fi).^2);

rs = 1 - (SSres/SStot);