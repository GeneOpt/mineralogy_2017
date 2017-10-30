% this script is used to generate I and P plots for all glaciers but this
% time the glaciers will be plotted together. 

clear all
close all
run loadSample_specs
run mineral_colors.m
    
labelR = {'1','2','4','5','6a','6b','7','8','9a','9b','10a','10b','11a',...
    '11b','14a','14b','16a','16b','17','18a','18b','19a','19b','20a','20b','21'}
pr = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0];

labelS = {'1','2','3','4','5','6','7','8','9','10','11','13',...
    '14','15','16','17','18','19','20','21','BH_1'}

s2rA = [1 NaN; 2 NaN; NaN NaN; 3 NaN; 4 NaN; 5 6; 7 NaN; 8 NaN; 9 10; 11 12;...
    13 14; NaN NaN; 15 16; 16 NaN; 17 18; 19 NaN; 20 21; 22 23; 24 25; 26 NaN]
% get all the folder names (each folder is a rock sample)

folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
[nmsR] = dir([folderR '\*.mat']);
matNmR = {nmsR.name}

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nmsS] = dir([folderS '\*.mat']);
matNmS = {nmsS.name}

dummyV = load([folderR 'GL01RS01.mat']);
ID = dummyV.isleD;
fID = fields(ID)

colR = [0.3 0.3 0.3;1 0 0];
colS = [0.7 0.7 0.7;1 0 1];
sS = {'-','--'}
sR = {'-','--'}


for M = 1:length(fID)
    close all
    f1 = figure
%     fn = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\islands\sizeDist\bins10\' matNmS{j}(1:5)]
%     mkdir(fn)
    
    for i = 1:length(matNmS)
        mns = matNmS{i};
        varS = load([folderS mns]);
        islS = varS.isleD;
        ptcS = varS.ptclD;
        gS = varS.D_c;
        mnrl_S = varS.mnrlMtx_c;
        [mnrlS,minsN,minNFull] = abbvMins(mnrl_S,5);  
        islS_M = islS.(fID{M});
        cnt = 1;
        elS = [];

        for k = 1:size(mnrlS,1)
            a = mnrlS(k,:);
            aN = a(M)./sum(a(1:38));
            if aN>0
                elS(cnt) = k;
                cnt = cnt+1;
            end

        end        
        D_mS = gS(elS);
        [binC,binC_S,yBinI,mBin,stdBin,sDat] = bin_szHist(islS_M,islS_M);
        [binC,binC_S,yBinP,mBin,stdBin,sDat] = bin_szHist(D_mS,D_mS);

        subplot(2,2,2)
        hold on
        yI = yBinI/sum(yBinI);
        yIMs(i) = max(yI);
        p2(i) = plot(binC,yI,sS{stS(i)+1},'color',colS(mxS(i)+1,:),'linewidth',1);
%         text(binC(1),yBinI(1)/sum(yBinI),labelS{i})
        xlabel('Island diameter - sediment (\phi)')
        subplot(2,2,4)
        hold on
        yP = yBinP/sum(yBinP);
        yPMs(i) = max(yP);
        p4(i) = plot(binC,yP,sS{stS(i)+1},'color',colS(mxS(i)+1,:),'linewidth',1);
        xlabel('Grain diameter - sediment (\phi)')
%         text(binC(end),yBinP(end)/sum(yBinP),labelS{i})
%         keyboard
    end

    
    for i = 1:length(matNmR)
        
        mns = matNmR{i};
        varR = load([folderR mns]);
        islR = varR.isleD;
        ptcR = varR.ptclD;
        gR = varR.D_c;
        mnrl_R = varR.mnrlMtx_c;
        [mnrlR,minsN,minNFull] = abbvMins(mnrl_R,5);  
        islR_M = islR.(fID{M});
        cnt = 1;
        elR = [];

        for k = 1:size(mnrlR,1)
            a = mnrlR(k,:);
            aN = a(M)./sum(a(1:38));
            if aN>0
                elR(cnt) = k;
                cnt = cnt+1;
            end
        end
        
        D_mR = gR(elR);
        [binC,binC_R,yBinI,mBin,stdBin,sDat] = bin_szHist(islR_M,islR_M);
        [binC,binC_R,yBinP,mBin,stdBin,sDat] = bin_szHist(D_mR,D_mR);
       
        subplot(2,2,1)
        hold on
        yI = yBinI/sum(yBinI);
        yIMr(i) = max(yI);
        p1(i) = plot(binC,yI,sR{stR(i)+1},'color',colR(pr(i)+1,:),'linewidth',1);
        xlabel('Island diameter - rock (\phi)')
%         text(binC(1),yBinI(1)/sum(yBinI),labelR{i})
        subplot(2,2,3)
        hold on
        yP = yBinP/sum(yBinP);
        yPMr(i) = max(yP);
        p3(i) = plot(binC,yP,sR{stR(i)+1},'color',colR(pr(i)+1,:),'linewidth',1);
        if i == 18
            plot(binC,yP,sR{stR(i)+1},'color',[1 0.5 0.5],'linewidth',1);
        end
        xlabel('Grain diameter - rock (\phi)')
%         text(binC(end),yBinP(end)/sum(yBinP),labelR{i})
%         keyboard
    end
        
    for sp = 1:4
        subplot(2,2,sp)
        grid on
        ylabel('Number of particles (pdf)');
        set(gca,'fontsize',15);
%         ylim([0 0.4])
        xlim([-10.5 -5.5])
    end
    
    yMu = max([yIMs,yIMr])
    subplot(2,2,1)
    text(0.8,0.9,minsN{M},'units','normalized','fontsize',20)
    ylim([0 yMu])
    subplot(2,2,2)
    ylim([0 yMu])
    
    yMl = max([yPMs,yPMr])
    subplot(2,2,3)
    ylim([0 yMl])
    subplot(2,2,4)
    ylim([0 yMl])
    
    l1 = legend([p1(2),p1(4),p1(18),p1(1),p2(2),p2(3),p2(14),p2(1)],...
        {'MSNS','MSS','PNS','PS','MSNS','MSS','MXNS','MXS'},'position',[0.8827    0.6621    0.0732    0.3109],...
        'units','normalized')

    
    saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\islands\sizeDist\bins5\all_together_now\'  fID{M}])
    savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\islands\sizeDist\bins5\all_together_now\'  fID{M}])

end

    

