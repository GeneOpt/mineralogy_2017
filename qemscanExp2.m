%% here is wahat I am doing to get a first pass at the QEMSCAN data.

clear all
close all

[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',2,'A1:AQ21');
samples = labs(2:end,1)
minerals = labs(1,6:end-5)
minDat = data(:,3:end-5)
glType = data(:,end)


f1 = figure(1)

for i = 1:length(samples)
    
    hold on
    [mlc, mfc, sz, tp] = specsCompare(glType(i))
    h(i) = plot(minDat(i,:),'marker', tp, 'markersize', sz,...
        'markerfacecolor', mfc, 'markeredgecolor', mlc, 'linestyle', 'none')
%     text(1:length(minDat(1,:)),minDat(i,:),repmat(samples(i),1,length(minDat)), 'fontsize', 8)

end

lh = legend([h(1),h(2),h(7),h(13),h(end)], {'MSNS',...
    'MSS','MXS','MXNS','GL1 BH'},...
    'fontsize', 7, 'location', 'northeast')  
grid on

XX = 1:length(minerals)
set(gca, 'XTick', XX, 'XTickLabel', '', 'fontsize', 16)
h = get(gca,'xlabel');
set(h, 'Units', 'data')
pos = get(h, 'position');
yy = pos(2);
for i = 1:length(minerals)
    t = text(XX(i),yy+0.2,minerals(i))
    set(t,'rotation', 70, 'horizontalalignment', 'right')
    set(t,'verticalalignment',...
        'middle', 'fontsize',10)
end
xlim([0 length(minerals)+1])
% savePDFfunction(f1,'BMMP_MI8118-Apr16')


[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',3,'A1:AQ21');
minDatA = data(:,3:end-5)
%% here you can plot the area percent against the mass percent 

fs = 14

for j = 1:length(minerals)
    
    close all
    f1 = figure
    x = minDat(:,j)
    y = minDatA(:,j)
    for i = 1:length(samples)

        hold on
        [mlc, mfc, sz, tp] = specsCompare(glType(i))
        h(i) = plot(x(i),y(i),'marker', tp, 'markersize', sz,...
            'markerfacecolor', mfc, 'markeredgecolor', mlc, 'linestyle', 'none')

    end
    xlabel('Mass %', 'fontsize', fs)
    ylabel('Area %', 'fontsize', fs)
    text(0.1,0.95,minerals(j),'fontsize', fs, 'units', 'normalized')
    legend([h(1),h(2),h(end-4),h(13),h(end)], {'MS-NS','MS-S','MX-S','MX-NS','GL1 BH'}, 'location', 'southeast')
    grid on
    lrf = refline(1,0)
    lrf.Color = 'k'
    print(f1,['D:\Code\Summer_2013_data\EXP2_Figures\qemscan\areaVsMass\' minerals{j}], '-djpeg')

end

%% perhaps check out the particle size distributions?
close all
f = figure
[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',4,'A1:H141');

samps4psd = {'2','3','4','5','6','7','8','9','10','11','12',...
             '13','14','15','16','17','18','19','20'}

sizes = [64,32,24,16,8,4,2,1]
sizesM = sizes(1:end-1)-(sizes(1:end-1)-sizes(2:end))/2
sizesM = sizesM(end:-1:1)/1000
for i = 1:length(glType)-1
    
    eyes = (1:7)+(7*(i-1))
    psd = data(eyes,4)
    psd = psd(end:-1:1)
    
    if glType(i) == 1
        subplot(2,2,3)
        col = [0.4 0.4 0.4]
        xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
        text(0.1,0.9,'MS-NS','units','normalized')
    end
    if glType(i) == 2
        subplot(2,2,1)
        col = [1 0 0]
        xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
        text(0.1,0.9,'MS-S','units','normalized')
    end
    if glType(i) == 3
        subplot(2,2,2)
        col = [1 0 0]
        xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
        text(0.1,0.9,'MX-S','units','normalized')
    end
    if glType(i) == 4
        subplot(2,2,4)
        col = [0.4 0.4 0.4]
        xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
        text(0.1,0.9,'MX-NS','units','normalized')
    end
    hold on
    plot(log2(sizesM),psd,'-o', 'color',col)
    text(log2(sizesM(6)),psd(6),samps4psd(i), 'fontsize',12)
    grid on
    
end
savePDFfunction(f,'QEM_GSD_4plot')
sizesMR = [0.02	0.488	0.977	1.953	3.906	5.86	7.8125	11.719	15.625	23.4375	31.25	46.875	62.5	78.125	93.75	109.375	125	156.25]

psdMR = [0	0.328015	3.915651	10.169406	22.121996	33.058425	42.724045	58.239745	69.459494	83.366187	90.750298	97.056098	99.093525	99.777382	99.98164	100	100	100]




%% compare masterszier vs qemscan for GL2
close all

[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',4,'A1:H141');

samps4psd = {'2','3','4','5','6','7','8','9','10','11','12',...
             '13','14','15','16','17','18','19','20'}

sizes = [64,32,24,16,8,4,2,1]
dxQ = sizes(1:end-1)-sizes(2:end)
sizesM = sizes(1:end-1)-(dxQ/2)
% for i = 1:length(glType)
figure
for i = 1
    
    eyes = (1:7)+(7*(i-1))
    psd = data(eyes,4)

    plot(log2(sizesM/1000),psd,'-o', 'color',col)
%     text(log2(sizesM(6)),psd(6),samps4psd(i), 'fontsize',12)
    grid on
    
end

sizesMR = [0.02	0.488	0.977	1.953	3.906	5.86	7.8125	11.719	15.625	23.4375	31.25	46.875	62.5	78.125	93.75	109.375	125	156.25]

psdMR = [0	0.328015	3.915651	10.169406	22.121996	33.058425	42.724045	58.239745	69.459494	83.366187	90.750298	97.056098	99.093525	99.777382	99.98164	100	100	100]

hold on
% plot(log2(sizesMR/1000),psdMR,'r-o')
dx = sizesMR(2:end)-sizesMR(1:end-1)
smrM = sizesMR(1:end-1)+(dx/2)
GSD_mrV = (psdMR(2:end)-psdMR(1:end-1))./sizesMRmid
plot(log2(smrM/1000),GSD_mrV,'r-^')


% convert mastersizer volume to area

GSD_mrA = (3/2)*GSD_mrV./smrM
hold on
plot(log2(smrM/1000),GSD_mrA,'r-s')

% compute the mean area grain size

x = sizesM(end:-1:1)/1000
y = psd(end:-1:1)

nx = 200
xq = linspace(x(1),x(end),nx);

dx = xq(3)-xq(2);
yq = interp1(x,y,xq,'pchip');

hold on
plot(log2(xq),yq)

A = sum(yq)*dx
yqN = yq/A
sum(yqN)*dx

hold on
plot(log2(xq),yqN)

mean = sum(yq.*xq)*dx
l2m = log2(mean)
xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')



%% compare masterszier means
close all

[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',4,'A1:H141');

sizes = [64,32,24,16,8,4,2,1]
dxQ = sizes(1:end-1)-sizes(2:end)
sizesM = sizes(1:end-1)-(dxQ/2)
x = sizesM(end:-1:1)/1000
nx = 200
xq = linspace(x(1),x(end),nx);
dx = xq(3)-xq(2);

figure

for i = 1:length(glType)
    
    eyes = (1:7)+(7*(i-1))
    psd = data(eyes,4)


    y = psd(end:-1:1)
    yq = interp1(x,y,xq,'pchip');

    subplot(4,5,i)
%     plot(log2(xq),yq)
    A = sum(yq)*dx
    yqN = yq/A
    sum(yqN)*dx
    hold on
    plot(log2(xq),yqN)
    mean(i) = sum(yq.*xq)*dx
    l2m = log2(mean(i))
    hold on
    plot([l2m l2m], get(gca,'Ylim'),'k')
    text(0.1,0.9,samples(i),'units','normalized')
    grid on
    ylim([0 60])
    
end


f = figure
for i = 1:length(glType)-1
    hold on
    [mlc, mfc, sz, tp] = specsCompare(glType(i));
    h(i) = plot(glType(i),mean(i),'marker', tp, 'markersize', sz,...
            'markerfacecolor', mfc, 'markeredgecolor', mlc,...
            'linestyle', 'none');
    text(glType(i),mean(i),samples{i})

end
xlabel('$$\rm log_2(mean\,grain\,size\,\,\mu m)$$','interpreter','latex')
    
xlim([0 5])
grid on
glString = {'MSNS','MSS','MXS','MXNS'}
set(gca, 'XTick', 1:length(glString), 'XTickLabel', glString, 'fontsize', 10);
        xlim([0 length(glString)+1]);   

savePDFfunction(f,'meanGSDfromQEM')   

%% here you want to plot things again the XRD data

[XRDsamps] = xlsread('D:\Field_data\2013\Summer\Geochemistry\EXP2_results.xlsx',2,'B1:U1')
[XRDdata,XRDmins] = xlsread('D:\Field_data\2013\Summer\Geochemistry\EXP2_results.xlsx',2,'A24:U40')

QSsamps = samples(1:end-1)
QSdata = minDat(1:end-1,:)
QSmins = minerals

XRDsamps = num2str(XRDsamps(2:end))
XRDdata = XRDdata(:,2:end)'

QSdat2 = [QSdata(:,28)+QSdata(:,29),...
            QSdata(:,3)+QSdata(:,4)+QSdata(:,5)+QSdata(:,6),...
            QSdata(:,19)+QSdata(:,20),...
            QSdata(:,8)+QSdata(:,9)+QSdata(:,10)+QSdata(:,11),...
            QSdata(:,18)...
            QSdata(:,13)+QSdata(:,14),...
            QSdata(:,23),...
            QSdata(:,30),...
            QSdata(:,16)+QSdata(:,17)+QSdata(:,7),...
            QSdata(:,30),...
            QSdata(:,2),...
            QSdata(:,1),...
            QSdata(:,22),...
            QSdata(:,14),...
            zeros(length(QSsamps),1),...
            QSdata(:,15),...
            zeros(length(QSsamps),1)]

size(XRDdata)
size(QSdat2)

fs = 16

% for j = 1:length(XRDmins)
for j = 1
    
    close all
    f1 = figure
    x = XRDdata(:,j)*100
    y = QSdat2(:,j)
    p = polyfit(x,y,1)
    m = p(1)
    b = p(2)
    
    for i = 1:length(XRDdata)

        hold on
        [mlc, mfc, sz, tp] = specsCompare(glType(i))
        h(i) = plot(x(i),y(i),'marker', tp, 'markersize', sz,...
            'markerfacecolor', mfc, 'markeredgecolor', mlc, 'linestyle', 'none')

    end
    xlabel('Talc (mass % XRD)', 'fontsize', fs)
    ylabel('Mg silicates (mass % Qemscan)', 'fontsize', fs)
%     text(0.9,0.15, XRDmins(j),'fontsize', fs, 'units', 'normalized')
    legend([h(1),h(2),h(7),h(13)], {'MS-NS','MS-S','MX-S','MX-NS'}, 'location', 'northwest')
    grid on
    lrf = refline(1,0)
    lrf.Color = 'k'
    xlm = get(gca,'xlim')
    hold on
    plot(xlm,m.*xlm+b,'b--')
    txt1 = text(0.2,0.9,['y = ' num2str(m,2) 'x + ' num2str(b,2) '\newline r = ' num2str(corr(x,y),2)],'units','normalized')
    txt1.FontSize = 16
    txt1.Color = 'b'
    txt2 = text(0.95,0.9,'1:1','units','normalized','fontsize',16)
    set(gca,'fontsize',16)
    savePDFfunction(f1,'talc')
    
    perd = (x-y)./y*100
    m_diff = mean(perd)
    stdp = std(perd)
%     print(f1,['D:\Code\Summer_2013_data\EXP2_Figures\qemscan\XRDvsQS\' XRDmins{j}], '-djpeg')
end


%% okay, now this is where the analysis gets involved, because I am going to 
% plot the mineralogy by grian size. plot them BY SIZE, but have the choice
% to group the minerals by like minerlas, and plot 2 suplots

[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',7,'A1:AR81')
minDat = data(:,5:end)
minerals = labs(1,8:end-5)

groupTheMins = 1
if groupTheMins == 1
    
   minDat = [minDat(:,1:2),minDat(:,3)+minDat(:,4)+minDat(:,5)+minDat(:,6),...
       minDat(:,7),minDat(:,8)+minDat(:,9)+minDat(:,10)+minDat(:,11),...
       minDat(:,12:15),minDat(:,16)+minDat(:,17),minDat(:,18:19),...
       minDat(:,20)+minDat(:,21),minDat(:,22:32)]
    
   minerals = {'Qtz','Kspar','Plag','Musc','Biot','Kaol','Chlr',...
    'Mg Clay','Mg Sil','Ill_Smec','Calc','Dolom','Fe Ox','Pyrite','Gyp_An',...
    'Halite','Rutile','Tita','Laum','Cpx','Fe Amph','Epidt','Apatite','Tourm'}
   
end

sizes = [64,32,16,4,1]
dx = sizes(1:end-1)-sizes(2:end)
sizeM = sizes(1:end-1)-dx/2

labels = [2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 21]

for j = 1:length(minerals_grouped)
    close all
    f1 = figure
    for i = 1:length(glType)-1

        eyes = (1:4)+(4*(i-1))
        psd = minDat(eyes,j)
    %     psd = psd(end:-1:1)

        if glType(i) == 1
            subplot(2,1,1)
            col = [0.4 0.4 0.4]
            xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MS-NS','units','normalized')
            ylabel('Mass %')
        end
        if glType(i) == 2
            subplot(2,1,1)
            title(minerals{j}, 'position', [0.95 1.05],'units','normalized')
            col = [1 0 0]
            xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
            text(0.1,0.9,'MS','units','normalized')
        end
        if glType(i) == 3
            subplot(2,1,2)
            col = [1 0 0]
            xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MX-S','units','normalized')
        end
        if glType(i) == 4
            subplot(2,1,2)
            col = [0.4 0.4 0.4]
            xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
            text(0.1,0.9,'MX','units','normalized')
        end
        hold on
        plot(log2(sizeM),(psd),'-o', 'color',col)
        text(log2(sizeM(2))+0.25,(psd(2)),num2str(labels(i)), 'fontsize',8)
        grid on
        ylabel('Mass %')

    end

print(f1,['D:\Code\Summer_2013_data\EXP2_Figures\qemscan\MinbySize_MassPercent_Group\bySize\2pgm\normal\' minerals{j}], '-djpeg')
end


%% okay, now this is where the analysis gets involved, because I am going to 
% plot the mineralogy by grian size. plot them BY MINERALOGY

[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',8,'A1:AR81')
minDat = data(:,5:end-5)

[s1 s2] = sort(mean(minDat))
s2 = s2(end:-1:1)
minDatS = minDat(:,s2)

sizes = [64,32,16,4,1]
dx = sizes(1:end-1)-sizes(2:end)
sizeM = sizes(1:end-1)-dx/2

sizeStr = {'32','32-16','16-4','4-1'}

labels = [2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 21]

minerals = minerals(s2)

for j = 1:length(sizeM)
    close all
    f1 = figure
    for i = 1:length(glType)-1

        psd = minDatS([1:4:4*length(samples)]+(j-1),:)
    %     psd = psd(end:-1:1)

        if glType(i) == 1
            sp1 = subplot(2,1,1);
            col = [0.4 0.4 0.4]            
        end
        if glType(i) == 2
            sp1 = subplot(2,1,1);
            title(sizeStr{j}, 'position', [0.95 1.05],'units','normalized')
            col = [1 0 0]
%             xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MS-S','units','normalized')

        end
        if glType(i) == 3
            sp2 = subplot(2,1,2);
            col = [1 0 0]
%             xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MX-S','units','normalized')
        end
        if glType(i) == 4
            sp2 = subplot(2,1,2);
            col = [0.4 0.4 0.4]
%             xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MX-NS','units','normalized')
        end
        hold on
        hp(i) = plot((psd(i,:)'),'color',col)
%         text(log2(sizeM(2))+0.25,psd(2),num2str(labels(i)), 'fontsize',8)
        grid on

    end

    XX = 1:length(minerals)
    
    subplot(2,1,1)
    set(gca, 'XTick', XX, 'XTickLabel', '', 'fontsize', 16)
    h = get(gca,'xlabel');
    set(h, 'Units', 'data')
    pos = get(h, 'position');
    yy = pos(2);
    for i = 1:length(minerals)
        t = text(XX(i),yy+0.2,minerals(i))
        set(t,'rotation', 70, 'horizontalalignment', 'right')
        set(t,'verticalalignment',...
            'middle', 'fontsize',10)
    end
    xlim([0 length(minerals)+1]) 
    text(0.9,0.9,'MS','units','normalized', 'fontsize',14)
    legend([hp(1) hp(2)], {'surge', 'non-surge'}, 'location', 'east', 'fontsize',10)
    ylabel('Area %')
    
    subplot(2,1,2)
    set(gca, 'XTick', XX, 'XTickLabel', '', 'fontsize', 16)
    h = get(gca,'xlabel');
    set(h, 'Units', 'data')
    pos = get(h, 'position');
    yy = pos(2);
    for i = 1:length(minerals)
        t = text(XX(i),yy+0.2,minerals(i))
        set(t,'rotation', 70, 'horizontalalignment', 'right')
        set(t,'verticalalignment',...
            'middle', 'fontsize',10)
    end
    xlim([0 length(minerals)+1]) 
    text(0.9,0.9,'MX','units','normalized', 'fontsize',14)
    ylabel('Area %')
    
print(f1,['D:\Code\Summer_2013_data\EXP2_Figures\qemscan\MinbySize_AreaPercent_Group\' sizeStr{j}], '-djpeg')
savePDFfunction(f1,['MinbySize_AreaPercent_Group\' sizeStr{j}])
end


%% okay, now this is where the analysis gets involved, because I am going to 
% plot the mineralogy by grian size. plot them BY MINERALOGY for the bulk
% mineralogy

[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',2,'A1:AQ21');
minerals = labs(1,6:end-5)

[data,labs] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8118-Apr16 ARQS Datasheet SFU Crompton.xlsx',2,'A1:AQ21');
samples = labs(2:end,1)
% minerals = labs(1,6:end-5)
minDat = data(:,3:end-5)
glType = data(:,end)


% if you want to group the minerals by similar minerlas, then run this
% section
gpMins = 1
if gpMins == 1

QSdata = minDat

minDat = [QSdata(:,28)+QSdata(:,29),...
            QSdata(:,3)+QSdata(:,4)+QSdata(:,5)+QSdata(:,6),...
            QSdata(:,19),...
            QSdata(:,8)+QSdata(:,9)+QSdata(:,10)+QSdata(:,11),...
            QSdata(:,18)...
            QSdata(:,13)+QSdata(:,14),...
            QSdata(:,23),...
            QSdata(:,30),...
            QSdata(:,16)+QSdata(:,17)+QSdata(:,7),...
            QSdata(:,30),...
            QSdata(:,2),...
            QSdata(:,1),...
            QSdata(:,22),...
            QSdata(:,14),...
            zeros(length(samples),1),...
            QSdata(:,15),...
            zeros(length(samples),1)]
[XRDdata,XRDmins] = xlsread('D:\Field_data\2013\Summer\Geochemistry\EXP2_results.xlsx',2,'A24:U40')

[s1 s2] = sort(mean(minDat))
s2 = s2(end:-1:1)
minDatS = minDat(:,s2)
minerals = XRDmins(s2)

end



% [s1 s2] = sort(mean(minDat))
% s2 = s2(end:-1:1)
% minDatS = minDat(:,s2)

sizes = [64,32,16,4,1]
dx = sizes(1:end-1)-sizes(2:end)
sizeM = sizes(1:end-1)-dx/2

sizeStr = {'32','32-16','16-4','4-1'}

labels = [2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 21]

% minerals = minerals(s2)

close all
f1 = figure
for i = 1:length(glType)-1

    psd = minDatS(:,:)
%     psd = psd(end:-1:1)

    if glType(i) == 1
        sp1 = subplot(2,1,1);
        col = [0.4 0.4 0.4]            
    end
    if glType(i) == 2
        sp1 = subplot(2,1,1);
        title('bulk mineralogy', 'position', [0.95 1.05],'units','normalized')
        col = [1 0 0]
%             xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MS-S','units','normalized')

    end
    if glType(i) == 3
        sp2 = subplot(2,1,2);
        col = [1 0 0]
%             xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MX-S','units','normalized')
    end
    if glType(i) == 4
        sp2 = subplot(2,1,2);
        col = [0.4 0.4 0.4]
%             xlabel('$$\rm log_2(grain\,size\,\,\mu m)$$','interpreter','latex')
%             text(0.1,0.9,'MX-NS','units','normalized')
    end
    hold on
    hp(i) = plot((psd(i,:)'),'color',col)
%         text(log2(sizeM(2))+0.25,psd(2),num2str(labels(i)), 'fontsize',8)
    grid on

end

XX = 1:length(minerals)

subplot(2,1,1)
set(gca, 'XTick', XX, 'XTickLabel', '', 'fontsize', 16)
h = get(gca,'xlabel');
set(h, 'Units', 'data')
pos = get(h, 'position');
yy = pos(2);
for i = 1:length(minerals)
    t = text(XX(i),yy+0.2,minerals(i))
    set(t,'rotation', 70, 'horizontalalignment', 'right')
    set(t,'verticalalignment',...
        'middle', 'fontsize',10)
end
xlim([0 length(minerals)+1]) 
text(0.9,0.9,'MS','units','normalized', 'fontsize',14)
legend([hp(1) hp(2)], {'surge', 'non-surge'}, 'location', 'east', 'fontsize',10)
ylabel('(Mass %)')

subplot(2,1,2)
set(gca, 'XTick', XX, 'XTickLabel', '', 'fontsize', 16)
h = get(gca,'xlabel');
set(h, 'Units', 'data')
pos = get(h, 'position');
yy = pos(2);
for i = 1:length(minerals)
    t = text(XX(i),yy+0.2,minerals(i))
    set(t,'rotation', 70, 'horizontalalignment', 'right')
    set(t,'verticalalignment',...
        'middle', 'fontsize',10)
end
xlim([0 length(minerals)+1]) 
text(0.9,0.9,'MX','units','normalized', 'fontsize',14)
ylabel('(Mass %)')
    
% print(f1,['D:\Code\Summer_2013_data\EXP2_Figures\qemscan\MinbySize_AreaPercent_Group\' sizeStr{j}], '-djpeg')
savePDFfunction(f1,'minMassPercent_2plot_groupedMins')






