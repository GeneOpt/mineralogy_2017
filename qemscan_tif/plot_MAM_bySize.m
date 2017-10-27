run loadSample_specs

varS = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\GL 01.mat')
varR = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\GL01B.mat')
MAM_cS = varS.MAM_c;
MAM_cR = varR.MAM_c;
MAM_all_cS = varS.MAM_all_c;
MAM_all_cR = varR.MAM_all_c;

AbM = 1

if AbM == 1
    [MAM_cS,minsN,minNFull] = abbvMins(MAM_cS,4)
    [MAM_cR,minsN,minNFull] = abbvMins(MAM_cR,4)
    [MAM_all_cS] = abbvMins(MAM_all_cS,3)
    [MAM_all_cR] = abbvMins(MAM_all_cR,3)
end

%% make an correlation plot of the MAM_all matrix for sed
f1 = figure
mcs = log(MAM_all_cS./repmat(sum(MAM_all_cS,2),[1,size(MAM_cS,1)])*100)
imagesc(mcs)
set(gca,'xtick',1:length(minsN),'xticklabel',minsN)
rotateXLabels(gca,60)
set(gca,'ytick',1:length(minsN),'yticklabel',minsN)
title('Glacier 1 sediment')
caxis([-6 4])
colorbar
% saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\jpg\GL01\GL01_corrAll_sed'])
% savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\pdf\GL01\GL01_corrAll_sed'])
%% make an correlation plot of the MAM_all matrix for rock
f2 = figure
mcr = log(MAM_all_cR./repmat(sum(MAM_all_cR,2),[1,size(MAM_cR,1)])*100)
imagesc(mcr)
set(gca,'xtick',1:length(minsN),'xticklabel',minsN)
rotateXLabels(gca,60)
set(gca,'ytick',1:length(minsN),'yticklabel',minsN)
title('Glacier 1 BH')
caxis([-6 4])
colorbar
saveJPGfunction(f2,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\jpg\GL01\GL01_corrAll_BH'])
savePDFfunction(f2,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\pdf\GL01\GL01_corrAll_BH'])
%% make a correalation plot of the difference
f3 = figure
imagesc(mcr-mcs)
set(gca,'xtick',1:length(minsN),'xticklabel',minsN)
rotateXLabels(gca,60)
set(gca,'ytick',1:length(minsN),'yticklabel',minsN)
title('Glacier 1 BH - sediment')
caxis([-6 4])
colorbar
saveJPGfunction(f3,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\jpg\GL01\GL01_corrAll_diffBmS'])
savePDFfunction(f3,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\pdf\GL01\GL01_corrAll_diffBmS'])

[~,isr] = sort(mean(MAM_all_cS),'descend');

for P = 1:29
    f1 = figure
%     set(f1, 'Position', [-1279, 508, 1280, 907])
    
    for S = 1:4
        tArr = MAM_cS(P,:,S)
        tArrS = tArr(isr);
        numPt = sum(tArrS);
        hold on
        y = log10(tArrS/numPt*100);
        h(S) = plot(y,'k-^','linewidth',S)
        grid on
        ylabel(minNFull{P})
        set(gca,'xtick',1:length(minsN),'xticklabel',minsN(isr))
        rotateXLabels(gca,60)
        
    end
    

    for S = 1:4
        tArr = MAM_cR(P,:,S);
        tArrS = tArr(isr);
        numPt = sum(tArrS);
        hold on
        hR(S) = plot(log10(tArrS/numPt*100),'b-^','linewidth',S)
    end
    legend([h(1) h(2) h(3) h(4) hR(1)],{'Seds - \phi>9.5','8.5<\phi<=9.5','7.5<\phi<=8.5','\phi<7.5','BH - \phi>9.5'})
    title('GL 01, MX-S')
    set(gca,'fontsize',19)
    
    if AbM == 1
        saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\jpg\GL01\BHSed\GL01_' mins{P} '_BHSed'])
        savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\mins_abbv\pdf\GL01\BHSed\GL01_' mins{P} '_BHSed'])
    end
    if AbM ~= 1
        saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\all_mins\jpg\GL01\BHSed\GL01_' mins{P} '_BHSed'])
        savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\all_mins\pdf\GL01\BHSed\GL01_' mins{P} '_BHSed'])
    end

end