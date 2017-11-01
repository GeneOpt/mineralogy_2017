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