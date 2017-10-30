function [m, b, yy,zx,zy] = linReg(x,y)

m = (sum((x-mean(x)).*(y-mean(y))))/(sum((x-mean(x)).^2));
b = mean(y) - m*mean(x);

yy = m*x+b;

zx = (x-mean(x))/std(x);
zy = (y-mean(y))/std(y);

R = sum(zx.*zy)/(length(x)-1) ; %sample correlation coefficient
rs = R*R     ;                  %coefficient of determination
sxy2 = var(y)*(1-rs);

sxy = (sum((x-mean(x)).*(y-mean(y)))/length(x)); %sample covariance
r = sxy/(std(x,1)*std(y,1))  %sample correlation coeffiecient
rs2 = r^2;

sydotx = sqrt(sum((y-yy).^2)/length(y)); % standard error of estimate
sydotxsq = sydotx^2;
sydotxsq2 = var(y,1)*(1-r^2);  % standard error of estimate square


rss = sum((y-yy).^2)
ev = sum((yy-mean(y)).^2)
rs2 = ev/sum((y-mean(y)).^2)  %r square = coefficient of determination
                              % is what fraction of total variation can be
                              % explaines by regression
                          