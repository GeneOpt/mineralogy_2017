function [yF,k,lbd] = fitW_J(y)

% close all
% plot(y,'-o')

kA = [0.1:0.01:5];
lbdA = [0.1:0.01:5];
x = 1:length(y);
        
for i =1:length(kA)
    k = kA(i);
    for j = 1:length(lbdA)
        lbd = lbdA(j);
        yF = (k/lbd)*((x/lbd).^(k-1)).*exp(-(x/lbd).^k);
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

k = kA(rw);
lbd = lbdA(cl);
yF = (k/lbd)*((x/lbd).^(k-1)).*exp(-(x/lbd).^k);
r = sqrt(sum((y-yF).^2));
% 
% hold on
% plot(x,yF,'r-o')
% text(0.9,0.9,num2str(r),'units','normalized')
% 
% figure
% [xx,yy] = meshgrid(kA,lbdA)
% surface(xx',yy',r2)
% xlabel('k')
% ylabel('lambda')

