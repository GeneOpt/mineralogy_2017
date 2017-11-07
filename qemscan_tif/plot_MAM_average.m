clear all
close all
run loadSample_specs
run mineral_colors.m

%% sample specs for rock

% samplesR = {'GL01_1','GL02_2','GL04_1','GL05_1','GL06_2','GL06_4','GL07_1',...
%     'GL08_2','GL09_A','GL09_B','GL10_1','GL10_2','GL11_1','GL11_2','GL14_2',...
%     'GL14_3','GL16_1','GL16_2','GL17_1','GL18_1','GL18_2','GL19_1','GL19_2',...
%     'GL20_1','GL20_2','GL21_1'}; the entire suite
samplesR = {'GL01_1','GL02_2','GL04_1','GL05_1','GL06_2','GL06_4','GL07_1',...
    'GL08_2','GL09_A','GL09_B','GL10_1','GL10_2','GL11_1','GL11_2','GL14_2',...
    'GL14_3','GL16_1','GL16_2','GL17_1','GL18_1','GL18_2','GL19_1','GL19_2',...
    'GL20_1','GL20_2','GL21_1'};
labelR = {'1','2','4','5','6a','6b','7','8','9a','9b','10a','10b','11a',...
    '11b','14a','14b','16a','16b','17','18a','18b','19a','19b','20a','20b','21'}
xValR = [16 1 2 6 7 7 8 17 3 3 4 4 9 9 18 18 12 12 13 19 19 20 20 14 14 15]

for i = 1:length(samplesR)
    keepLabR(i) = find(strcmp(varsR(:,8),samplesR(i)))
end

stR = strcmp(varsR(keepLabR,2),'S')
nstR = strcmp(varsR(keepLabR,2),'NS')
gtR = varsR(keepLabR,2)
prR = strcmp(varsR(keepLabR,3),'P')
msR = strcmp(varsR(keepLabR,3),'MS')
groupR = datR(keepLabR-1,7)


%% sample specs for sediment
% samplesS = {'GL 01','GL 02','GL 03','GL 04','GL 05','GL 06','GL 07','GL 08',...
%     'GL 09','GL 10','GL 11','GL 13','GL 14','GL 15','GL 16','GL 17','GL 18',...
%     'GL 19','GL 20','GL 21','GL1_1'}; the entire suite
samplesS = {'GL 01','GL 02','GL 03','GL 04','GL 05','GL 06','GL 07','GL 08',...
    'GL 09','GL 10','GL 11','GL 13','GL 14','GL 15','GL 16','GL 17','GL 18',...
    'GL 19','GL 20','GL 21','GL1_1'};
labelS = {'1','2','3','4','5','6','7','8','9','10','11','13',...
    '14','15','16','17','18','19','20','21','BH_1'}
xValS = [16 1 5 2 6 7 8 17 3 4 9 10 18 11 12 13 19 20 14 15 16]   


for i = 1:length(samplesS)
    keepLabS(i) = find(strcmp(varsS(:,6),samplesS(i)))
end

stS = strcmp(varsS(keepLabS,2),'S')
nstS = strcmp(varsS(keepLabS,2),'NS')
gtS = varsS(keepLabS,2)
mxS = strcmp(varsS(keepLabS,3),'MX')
msS = strcmp(varsS(keepLabS,3),'MS')
groupS = datS(keepLabS-1,5)
labelS = labelS(keepLabS-1)
%%

folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
[nmsR] = dir([folderR '\*.mat']);
matNmR = {nmsR.name}

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nmsS] = dir([folderS '\*.mat']);
matNmS = {nmsS.name}


% min1 = find(strcmp(minsN,'Epd Zo'))
% min2 = find(strcmp(minsN,'Ab'))
minL = minsN(2:end);

%%
% for min1 = 2:length(minsN)
for min1 = 8:11
    
    for Hs = 1:length(xValS)    % for the sediment

        mns = matNmS{Hs};
        varS = load([folderS mns]);

    %     MAM_cS = varS.MAM_c;
        MAM_all_cS = varS.MAM_all_c;

        tArrS = MAM_all_cS(min1,2:end);
        tArrSP(Hs,:) = tArrS/sum(tArrS);
% return
    end
    %%

    f1 = figure
    fs = 7
    lw = 2

    subplot(2,2,3)
    yL(3) = ylabel('f.n. (MS - sediment)')
    hold on
    mtA = nanmean(tArrSP(msS,:));
    stdA = std(tArrSP(msS,:));
    [y,isr] = sort(mtA,'descend');
    stdY = stdA(isr);
    % y = log10(y)
    plot(1:length(y),y,'-','color',[0.5 0.5 0.5],'linewidth',lw);

    set(gca,'xtick',1:length(y),'xticklabel',minL(isr))
    rotateXLabels(gca(), 60)
    for s = 1:length(stdY)
        plot([s s],[y(s)+stdY(s) y(s)-stdY(s)],'k','linewidth',lw)
    end
    grid on
    set(gca,'fontsize',fs)

    subplot(2,2,4)
    yL(4) = ylabel('f.n. (MX - sediment)')
    hold on
    mtA = nanmean(tArrSP(mxS,:));
    stdA = std(tArrSP(mxS,:));
    [y,isr] = sort(mtA,'descend');
    stdY = stdA(isr);
    % y = log10(y)
    plot(1:length(y),y,'m-','linewidth',lw);
    % minL = minsN(2:end);
    set(gca,'xtick',1:length(y),'xticklabel',minL(isr))
    rotateXLabels(gca(), 60)
    for s = 1:length(stdY)
        plot([s s],[y(s)+stdY(s) y(s)-stdY(s)],'k','linewidth',lw)
    end
    grid on
    set(gca,'fontsize',fs)

    %% rock
    for Hs = 1:length(xValR)    % for the sediment

        mns = matNmR{Hs};
        varS = load([folderR mns]);

    %     MAM_cS = varS.MAM_c;
        MAM_all_cS = varS.MAM_all_c;

        tArrS = MAM_all_cS(min1,2:end);
        tArrSP(Hs,:) = tArrS/sum(tArrS);

    end
    %%

    subplot(2,2,1)
    yL(1) = ylabel('f.n. (MS - Rock)')
    hold on
    mtA = nanmean(tArrSP(msR,:));
    stdA = std(tArrSP(msR,:));
    [y,isr] = sort(mtA,'descend');
    stdY = stdA(isr);
    % y = log10(y)
    plot(1:length(y),y,'k-','linewidth',lw);
    % minL = minsN(2:end);
    set(gca,'xtick',1:length(y),'xticklabel',minL(isr))
    rotateXLabels(gca(), 60)
    for s = 1:length(stdY)
        plot([s s],[y(s)+stdY(s) y(s)-stdY(s)],'k','linewidth',lw)
    end
    grid on
    set(gca,'fontsize',fs)


    subplot(2,2,2)
    yL(2) = ylabel('f.n. (MX - Rock)')
    hold on
    mtA = nanmean(tArrSP(prR,:));
    stdA = std(tArrSP(prR,:));
    [y,isr] = sort(mtA,'descend');
    stdY = stdA(isr);
    % y = log10(y)
    plot(1:length(y),y,'r-','linewidth',lw);
    % minL = minsN(2:end);
    set(gca,'xtick',1:length(y),'xticklabel',minL(isr))
    rotateXLabels(gca(), 60)
    for s = 1:length(stdY)
        plot([s s],[y(s)+stdY(s) y(s)-stdY(s)],'k','linewidth',lw)
    end
    grid on
    set(gca,'fontsize',fs)

    for sp = 1:4
        subplot(2,2,sp)
        yL(sp).FontSize = 14
        ylim([0 0.4])
    end

    subplot(2,2,1)
    t1 = title(minNFull{min1},'fontsize',14)
    t1.Position = [34.3830    0.4241         0]

    saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\subplot4\' minL{min1-1}])
    savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\subplot4\' minL{min1-1}])
    close all
end
