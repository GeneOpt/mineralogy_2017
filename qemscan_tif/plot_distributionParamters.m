clear all
close all
run loadSample_specs
run mineral_colors.m

%% for the sediment

folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
dummyV = load([folderR 'GL01RS01.mat']);
ID = dummyV.isleD;
fID = fields(ID)

varLabs = {'m (fractal)','b (fractal)','r (fractal)',...
    'k (Weibull)', '\lambda (Weibull)','r (Weibull)',...
    '\mu (log normal)', '\sigma (log normal)','r (log normal)'}

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\fitGSD\islands\sed\'

[nms] = dir([folderS '\*.mat']);
matNmS = {nms.name}

load([folderS matNmS{1}])

[fnms] = dir([folderS]);
fNmS = {fnms.name}
fNmS = fNmS(3:23)

p1 = find(strcmp(varLabs,'r (fractal)'))
p2 = find(strcmp(varLabs,'r (Weibull)'))

M = find(strcmp(fID,'Ab'))

for i = 1:length(fNmS)
    x1(i) = paramM(M,p1,i);
    x2(i) = paramM(M,p2,i);
end

%%
xlm = 0.2
ylm = 0.2

colA = [0.7 0.7 0.7;1 0.7 0.7]
mt = {'^','o'}

f1 = figure
subplot(1,2,2)
for i = 1:length(fNmS)

    hold on
    plot(x1(i),x2(i),['k' mt{mxS(i)+1}],'markerfacecolor', colA(stS(i)+1,:))
    text(x1(i),x2(i),labelS{i})

end
title(['Sediment'])
grid on 
xlim([0 xlm])
ylim([0 ylm])
xlabel(varLabs{p1})
ylabel(varLabs{p2})
%% for the rock 

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\fitGSD\islands\rock\'

[nms] = dir([folderS '\*.mat']);
matNmS = {nms.name}

load([folderS matNmS{1}])

[fnms] = dir([folderS]);
fNmS = {fnms.name}
fNmS = fNmS(3:28)


for i = 1:length(fNmS)
    x1(i) = paramM(M,p1,i);
    x2(i) = paramM(M,p2,i);
end

%%
colA = [0.1 0.1 0.1;1 0 0]
mt = {'o','^'}

subplot(1,2,1)
for i = 1:length(fNmS)

    hold on
    plot(x1(i),x2(i),['k' mt{msR(i)+1}],'markerfacecolor', colA(stR(i)+1,:))
    text(x1(i),x2(i),labelR{i})

end
title(['Rock - ' fID{M}])
grid on 
xlim([0 xlm])
ylim([0 ylm])
xlabel(varLabs{p1})
ylabel(varLabs{p2})
