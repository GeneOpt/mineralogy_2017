load([folderR matNmFlR{k}]);
for m = 2:length(mins)
q = clustM_c.(mins{m});
qs = sqrt(q);
el0 = find(q(:,1)==0);
qs(el0,:) = [];
close all
figure
plot(qs(:,1),qs(:,2),'.')

bins = 5
[binC,yBin] = bin_szMean(qs(:,1),qs(:,2),bins)
hold on
plot(binC,yBin,'r.-')

title(mins(m))
refline(1,0)
keyboard
end

%% attempts at homebrew fourier
figure
subplot(2,1,1)
xx_1 = xx - xx(1);
plot(xx_1,yy,'k-o');

lx = length(xx)
jf = 0:lx-1; 
kArr = [1:20]
A_k = []
B_k = []

for kI = 1:length(kArr)
    k = kArr(kI)
    B_k(kI) = (2/k)*sum(yy.*sin(2*pi*jf*k/lx))
    A_k(kI) = (2/k)*sum(yy.*cos(2*pi*jf*k/lx))
end
A_k = [,A_k]
B_k = [B_k]
Amp_k = sqrt((A_k.^2)+(B_k.^2))
Amp_k = Amp_k/max(Amp_k)


A_k = A_k/max(A_k)
B_k = B_k/max(B_k)

subplot(2,1,2)
plot([kArr],Amp_k,'b+-')

th = linspace(0,2*pi,50)
T = th(2)-th(1)
Fs = 1/T
Y = zeros(1,length(th))
for N = 1:length(kArr)
    Y = (A_k(N)*cos(N*th) + B_k(N)*sin(N*th))+Y
end
subplot(2,1,1)
hold on
plot(th,Y+mean(R_j),'m')


N = length(xx_1) 
for p = 1:N/2+1
    A1(p) = 0;
    B1(p) = 0;
    for n = 1:N
        A1(p) = A1(p)+2/N*cos(2*pi*(p-1)*n/N);
        B1(p) = B1(p)+2/N*sin(2*pi*(p-1)*n/N);
    end
end
A1(N/2+1) = A1(N/2+1)/2;

C1 = sqrt((A1.^2)+(B1.^2))
pmax = 10
for n = 1:N
    ynew(n) = A1(1)/2;
    for p = 2:pmax
        ynew(n) = ynew(n)+A(p)*cos(2*pi*(p-1)*n/N)+B(p)*sin(2*pi*(p-1)*n/N);
    end
end
hold on
plot(xx_1,ynew,'r')

return