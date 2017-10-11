clear all

run('GL1_1214.m')
run('GL1_1217.m')
run('GL1_1225.m')

minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi'])


altMtx = [alt14;alt17;alt25]
minCount = sum(1 - (altMtx == 0))
altAll = (sum(altMtx)./minCount)

modalComp = [modalComp14;modalComp17;modalComp25]
avemodalComp = sum(modalComp)/3
modalComp = [modalComp;avemodalComp]


%% for plotting
close all
figure('units','normalized','outerposition',[0 0 1 1])

for i = 1:4
    if i == 1
        col = 'r'
        mark = 's'
        ms = 8
    end
    if i == 2
        col = 'm'
        mark = 's'
        ms = 8
    end
    if i == 3
        col = 'b'
        mark = 's'
        ms = 8
    end
    if i == 4
        col = 'k'
        mark = 'o'
        ms = 12
    end
    
    plot(1:length(minerals),modalComp(i,:), col, 'marker',...
        mark, 'markerfacecolor', col, 'markersize', ms, 'linestyle', 'none')
    hold on
end

ylim([-2 40])
xlim([0 10])
set(gca, 'XTick', 1:length(minerals), 'XTickLabel', minerals, 'fontsize', 14)
xl = xlabel('Mineral')
yl = ylabel('Modal percent')

set(xl, 'fontsize', 22)
set(yl, 'fontsize', 22)


%% for IUGS classification
IUGS = modalComp(4,[1,2,3])
totalSil = sum(IUGS)
IUGS = IUGS/totalSil

%% re-organize to plot with SSC
% TSModalComp = zeros(4,9)
TSModalComp(:,1)  = modalComp(:,6)
TSModalComp(:,2)  = modalComp(:,2)
TSModalComp(:,4)  = modalComp(:,7)
TSModalComp(:,6)  = modalComp(:,5)
TSModalComp(:,8)  = modalComp(:,9)
TSModalComp(:,11) = modalComp(:,3)
TSModalComp(:,13) = modalComp(:,1)
TSModalComp(:,14) = modalComp(:,4)


save('TSminComp.mat', 'TSModalComp')
minerals2 =     [{'Actinolite  '};{'Albite Low  '};...
{'Ankerite-Dol'};{'Biotite     '};{'Calcite     '};{'Chlinoclore '};...
{'Gypsum      '};{'Clinozoisite'};{'Illite/Mcsvt'};...
{'Laumontite  '};{'K-Feldspar  '};{'Montmorillon'};{'Quartz      '};...
{'Sphene      '}];

minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi'])




