close all
clear variables

M1 = 0.05
M2 = 0.6

C = 1-(M1+M2);

X1 = [0:0.1:2]
X2 = [1:0.1:2]
for i = 1:length(X2)
    
    frac1 = (X1.*(M1+M2+C))./(X1.*M1+ X2(i).*M2 + C)
    hold on
    title(num2str(X2(i)))
    plot(X1,frac1)
    grid on
    keyboard

end

