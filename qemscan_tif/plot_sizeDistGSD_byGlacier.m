clear all
close all
run loadSample_specs
run mineral_colors.m

fOut = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\fitGSD\'
% mkdir(fOut)

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

for M = 1
    f1 = figure
    
    for i = 1:length(SPIR)
        varR = load([folderR matNmR{i}]);
        mtx_M = varR.mnrlMtx_c;
        D = varR.D_c;
        [mtxM,minsN,minNFull] = abbvMins(mtx_M,5);
%         M_i = mtxM(:,M)./sum(mtxM(:,1:39),2);
        M_i = mtxM(:,M);
        islD = varR.isleD.(fID{M});
        
        subplot(6,4,SPIR(i))
        hold on
        
%         [binC,yBin] = bin_szHist_mnrlMtx(D,mtxM,M); % this is if normalizing by total of all minerals in a bin
%          y = yBin;
%         
%         [binC,binC_S,yBin,mBin,stdBin,sDat,sBin] = bin_szHist(D,M_i);
%         %this is if normalizing by mineral across all sizes
%         y = (sBin/sum(sBin))*100;

        [binC,binC_S,yBin,mBin,stdBin,sDat,sBin] = bin_szHist(D,D);
        %this is if using the island distribution plots
        y = ((yBin/sum(yBin))*100);



        p1(i) = plot(binC,y,'-','color',colR(prR(i)+1,:),'linewidth',1);
        text(binC(1),y(1),labelR{i})
        if i == 18
            p1(i) = plot(binC,y,'-','color',[1 0.5 0.5],'linewidth',1);
            text(binC(1),y(1),labelR{i})
        end      
        yMr(i)=max(y);
    end

    
    for i = 1:length(matNmS)
            

        varS = load([folderS matNmS{i}]);
        mtx_M = varS.mnrlMtx_c;
        D = varS.D_c;
        [mtxM,minsN,minNFull] = abbvMins(mtx_M,5);
%         M_i = mtxM(:,M)./sum(mtxM(:,1:39),2);
        M_i = mtxM(:,M);
        islD = varS.isleD.(fID{M});

  
        subplot(6,4,SPIS(i))
        if SPIS(i) == 1 || SPIS(i) == 9 || SPIS(i) == 18 || SPIS(i) == 22
            ylabel('Percent')
        end
        if SPIS(i) == 13 || SPIS(i) == 22 || SPIS(i) == 19 || SPIS(i) == 20
            xlabel('log[ Grain size (-\phi)]')
        end
        if SPIS(i) == 1 
            title(['MSNS - GSD'])
        end
        if SPIS(i) == 2 
            title('MSS')
        end
        if SPIS(i) == 3 
            title('MXNS')
        end
        if SPIS(i) == 4 
            title('MXS')
        end
        
        hold on
        
%         
%         [binC,yBin] = bin_szHist_mnrlMtx(D,mtxM,M); % this is if normalizing by total of all minerals in a bin
%         y = yBin;
        
%         [binC,binC_S,yBin,mBin,stdBin,sDat,sBin] = bin_szHist(D,M_i);
%         this is if normalizing by mineral across all sizes
%         y = (sBin/sum(sBin))*100;

        [binC,binC_S,yBin,mBin,stdBin,sDat,sBin] = bin_szHist(D,D);
        %this is if using the island distribution plots
        y = ((yBin/sum(yBin))*100);
        
        p2(i) = plot(binC,y,'--','color',colS(mxS(i)+1,:),'linewidth',1)    
        text(binC(1),y(1),labelS{i})
        grid on
        yMs(i) = max(y);
    end

    yMax = max([yMr, yMs]);
    
    for sp = [1:16,18:20,22]
        subplot(6,4,sp)
        ylim([0 yMax])
        xlim([-10.25 -5.75])
    end
    l1 = legend([p1(2) p1(1) p2(2) p2(1)], {'MS rock','P rock','MS sed','MX sed'},'position',[0.5648    0.1300    0.0446    0.0691])

    saveJPGfunction(f1,[fOut '_GSD'])
    savePDFfunction(f1,[fOut '_GSD'])

    close all

end
    