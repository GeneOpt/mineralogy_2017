clear all
close all
clc

k = 0.5
c1 = [1 1 0 0 0;...
      1 0.5 0.5 0 0;...
      1 0.5 0.25 0.25 0;...
      1.5 0.5 0 0 0;...
      1 1 1 1 1;...
      1 1 1 1 0;
      10 1 1 1 1;
      10 10 1 0 0;...
      1 1 1 0 0 ;...
      1 1 1 0.5 0.5];
  
c = c1./repmat(sum(c1,2),[1,size(c1,2)])

for i=1:size(c,1)
    
    a = sort(c(i,:),'descend')
   
    H(i) = a(1)-(a(2)*a(3))-(a(3)*a(4))-(a(4)*a(5));
    tic
    t(1) = a(1)
    for j = 1:size(c,2)-2
        t(j+1) = t(j)-(a(j+1)*a(j+2))
    end
    T(i) = t(end)
    t1 = toc
    S(i) = skewness(a)
    
    
    A(i) = a(1)+sum([a(2:end-1).*a(3:end)]*-1)
    t2 = toc
    t1
    t2-t1
  
    
end

f1 = figure
p1 = plot(H,S,'bo','markerfacecolor','b','markersize',10)
for i = 1:length(H)
    text(H(i),S(i),num2str(i))
end
grid on
ylabel('skewness')
xlabel('H = f_1 - ... - f_{n-1}*f_{n1}')
[[1:size(c1)]',c]
set(gca,'fontsize',26)
% savePDFfunction(f1,'D:\Documents\Writing\Thesis\Current_use\figures\rockFrags\homogeneity')
