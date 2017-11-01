function [yF,mu, s] = fintLN_J(y)

% close all
% plot(y,'-o')

sA = [0.1:0.01:5];
muA = [0.1:0.01:5];
x = 1:length(y);
        
for i =1:length(sA)
    s = sA(i);
    for j = 1:length(muA)
        mu = muA(j);
        yF = 1/(sqrt(2*pi*s^2))*(exp(-(x-mu).^2/(2*s^2)));
        r2(i,j) = sqrt(sum((y-yF).^2));
%         pause(0.1)
%         hold on
%         plot(x,yF,'r-o')
    end
end

elM = find(r2==min(min(r2)));
[rc] = ind2rc(size(r2,1),elM);
rw = int32(rc(1));
cl = int32(rc(2));

s = sA(rw);
mu = muA(cl);
yF = 1/(sqrt(2*pi*s^2))*(exp(-(x-mu).^2/(2*s^2)));
% r = rsquare(yF,y)
% 
% hold on
% plot(x,yF,'r-o')
% text(0.9,0.9,num2str(r),'units','normalized')
% 
% figure
% surface(r2)
% return
