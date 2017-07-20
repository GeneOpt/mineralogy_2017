close all
clear variables

run mineral_colors
mins = mins(2:end);
mins_p = mins

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\neighbourPixels\';
BHM = xlsread([folder 'GL1_13BH50.xlsx'],3,'B2:AE31');
sedM = xlsread([folder 'GL 01.xlsx'],3,'B2:AE31');
rockM = xlsread([folder 'GL1_12_25.xlsx'],3,'B2:AE31');

ind_min = 'Ill Smec'
ind_minIND = (find(strcmp(mins,ind_min)))
mins_p(ind_minIND) = []

f1 = figure
hold on
xArr = 1:length(mins)
xArr(ind_minIND) = []
yBH = BHM(xArr,ind_minIND);
ySed = sedM(xArr,ind_minIND);
yRock = rockM(xArr,ind_minIND);

[srt,isr] = sort(yBH,'descend');

p1 = plot(xArr,log(yBH(isr)),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','k','markersize',10)
p2 = plot(xArr,log(ySed(isr)),'s','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k','markersize',10)
p3 = plot(xArr,log(yRock(isr)),'^','MarkerFaceColor',[0.6 0.6 0],'MarkerEdgeColor','k','markersize',10)

ylabel(['log percent of mineral neighbouring ' ind_min])
grid on
legend([p1 p2 p3],{'borehole','sediment','rock'},'location','northeast')
set(gca,'xtick',1:length(mins_p),'xticklabel',mins_p(isr))
rotateXLabels(gca,60)
set(gca,'fontsize',16)
savePDFfunction(f1, ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\percNeighb_' ind_min '_log'])