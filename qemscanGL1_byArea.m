

initz = 1
if initz == 1
    minPBulk = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\qemscan12Jul13_byArea.xlsx',2,'H2:AL2');
    minPSize = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\qemscan12Jul13_byArea.xlsx',3,'H2:AL6');
    minPSize = flipud(minPSize)/100;
    minPSiGr = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\qemscan12Jul13_byArea.xlsx',4,'H2:AL6');
    minPSiGr = flipud(minPSiGr)/100;
    [blank, names] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\qemscan12Jul13_byArea.xlsx',3,'H1:AL1');
end

sizes = [0.244,1,4,16,32,64];
% sizes = fliplr(sizes)
sMid = (sizes(2:end)+sizes(1:end-1))/2;


close all

load('GL1Mean.mat')
xgl1 = log2(x2d/1000);
ygl1 = meanGSD;
    
x = log2(sizes/1000);
% for i = 1:length(minPSize(1,:))
for i = 1:5
    
    figure
    M = i;
    y = minPSize(:,i);
        for j = 1:length(sizes)-1
            xf = [x(j) x(j+1) x(j+1) x(j)];
            yf = [0 0 y(j) y(j)];
            hold on
            fill(xf,yf,'b')
            hold on
            plot(xgl1,ygl1,'r','linewidth',2)
        end
%         hold on
%         bar(log2(sMid/1000),y,'m')

    title(names(i))
    xlim([-12 -2])

end

figure
y = sum(minPSize,2);
    for j = 1:length(sizes)-1
        xf = [x(j) x(j+1) x(j+1) x(j)];
        yf = [0 0 y(j) y(j)];
        hold on
        l1 = fill(xf,yf,'b');
        hold on
        l2 = plot(xgl1,ygl1,'r','linewidth',2);
    end
xlim([-12 -2])
title('All minerals')
xlabel('log2(size/mm)')
ylabel('pdf')
legend([l1 l2], {'qemscan','average of GL1 Jul 2013'},'location', 'northwest')

