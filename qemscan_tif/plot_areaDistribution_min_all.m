clear all
close all
run loadSample_specs
run mineral_colors.m

fOut = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\areaDistributions\bins5\allTogetherNow\'
mkdir(fOut)

folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
[nmsR] = dir([folderR '\*.mat']);
matNmR = {nmsR.name}

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nmsS] = dir([folderS '\*.mat']);
matNmS = {nmsS.name}

colR = [0.3 0.3 0.3;1 0 0];
colS = [0.7 0.7 0.7;1 0 1];
sS = {'-','--'}
sR = {'-','--'}

for M = 1:44
    f1 = figure
    subplot(1,2,1)
    hold on
    for i = 1:length(matNmR)
        varR = load([folderR matNmR{i}]);
        mtx_M = varR.mnrlMtx_c;
        D = varR.D_c;
        [mtxM,minsN,minNFull] = abbvMins(mtx_M,5);
%         M_i = mtxM(:,M)./sum(mtxM(:,1:39),2);
        M_i = mtxM(:,M);
        [binC,binC_S,yBin,mBin,stdBin,sDat,sBin] = bin_szHist(D,M_i);
        y = (sBin/sum(sBin));
        p1(i) = plot(binC,y,sR{stR(i)+1},'color',colR(prR(i)+1,:),'linewidth',1);
        if i == 18
            p1(i) = plot(binC,y,sR{stR(i)+1},'color',[1 0.5 0.5],'linewidth',1);
        end      
    end
    title('Rock')
%     text(0.1,0.9,minsN{M})
    grid on
    xlabel('Grain size (-\phi)')
    ylabel(['Area of ' minsN{M} ' (pdf)'])
    set(gca,'fontsize',18)
    xlim([-10.5 -5.5]);
    yl1 = get(gca,'ylim');
    
    subplot(1,2,2)
    hold on
    for i = 1:length(matNmS)
        varS = load([folderS matNmS{i}]);
        mtx_M = varS.mnrlMtx_c;
        D = varS.D_c;
        [mtxM,minsN,minNFull] = abbvMins(mtx_M,5);
%         M_i = mtxM(:,M)./sum(mtxM(:,1:39),2);
        M_i = mtxM(:,M);
        [binC,binC_S,yBin,mBin,stdBin,sDat,sBin] = bin_szHist(D,M_i);
        y = (sBin/sum(sBin));
        p2(i) = plot(binC,y,sS{stS(i)+1},'color',colS(mxS(i)+1,:),'linewidth',1)    
    end
    title('Sediment')
    xlabel('Grain size (-\phi)')
    grid on
    xlim([-10.5 -5.5])
    set(gca,'fontsize',18)
    l1 = legend([p1(2),p1(4),p1(18),p1(1),p2(2),p2(3),p2(14),p2(1)],...
        {'MSNS','MSS','PNS','PS','MSNS','MSS','MXNS','MXS'},'position',[0.1623    0.6033    0.0805    0.3724],...
        'units','normalized');
    yl2 = get(gca,'ylim');
    yMaxP = max([yl1(2) yl2(2)]);
    ymP(M) = yMaxP;
    
    for sp = 1:2
        subplot(1,2,sp)
        ylim([0 yMaxP])
    end
    
    return
    saveJPGfunction(f1,[fOut minsN{M}])
    savePDFfunction(f1,[fOut minsN{M}])
    
    close all
    
end
    