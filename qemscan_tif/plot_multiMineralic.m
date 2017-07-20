% in this script I am taking all of the matrices generated from
% "find_isolatedGrains", which is a matrix (mnrlMtx) of the number of pixels of each
% mineral (column) for each grain (row). The scale for a grain is 25/34
% microns per pixel
clc
close all
clear variables

run mineral_colors
[datR,varsR] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',1,'B1:I27');
[datS,varsS] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',2,'B1:I27');
RSArr = [1,2,4,5,6,7,8,9,9,10,10,11,11,13,13,15,15,16,17,17,18,18,19,19,20];
SPIS =  [4,1,2,5,6,10,14,8,9,13,18,22,12,3,7,11,16,20,15,19,4];
SPIR = [4,1,5,6,10,10,14,8,9,9,13,13,18,18,12,12,7,7,11,16,16,20,20,15,15,19];

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

%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder

sect = 0
if sect == 1

% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\seds\';
nmsFol = dir(folder)
matNmFl = {nmsFol.name}
matNmFl = matNmFl(3:end)

% grab the folder and look inside
for F = 1:length(matNmFl)
    inFold = [folder matNmFl{F}]
    [nms] = dir([inFold '\*.mat'])
    matNm = {nms.name}
    T_Arr = []
    S_Arr = []
    sz_Arr = []
    
    close all
    for N = 1:length(matNm)

        fname = [inFold '\' matNm{N}]

        load(fname)
        T = zeros(1,size(mnrlMtx,1));
        sz = zeros(1,size(mnrlMtx,1));
        if N == 1
            MnrlMtx = mnrlMtx;
        elseif N >1
            MnrlMtx = [MnrlMtx;mnrlMtx];
        end
        

        for i = 1:size(mnrlMtx,1)
            grain1 = sort(mnrlMtx(i,:),'descend');
            grain = grain1/sum(grain1);
            t = zeros(1,size(mnrlMtx,2)-1);
            t(1) = grain(1);
            for j = 2:size(mnrlMtx,2)-1
                t(j) = t(j-1)-(grain(j)*grain(j+1));
            end

            T(i) = t(end);
            sz(i) = sum(grain1); % size in pixels
            S(i) = skewness(grain);
        end

        T_Arr = [T_Arr,T];
        S_Arr = [S_Arr,S];
        sz_Arr = [sz_Arr,sz];

        f1 = figure(1)
        hold on
        plot(sz,T,'.')
        xlabel('grain size (pixels)')
        ylabel('homogeneity')

        % hold on
        % n = [1:500];
        % y = (n-2)./n
        % plot(n,y,'r-+')
        % 
        % y2 = (1-(1./n)+n)./(2*n-1)
        % plot(n,(1./y2)-1,'m-')


%         f2 = figure(2)
%         hold on
%         plot(sz,S,'.')
%         xlabel('grain size (pixels)')
%         ylabel('skewness')

    end
    
    mn = matNm{1}
    mn5 = mn(1:5)
    save([inFold '.mat'],...
        'T_Arr','S_Arr','sz_Arr','MnrlMtx')
    
end

end

%% loop through all the files and put the arrays into a matrix SED 
folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\seds\';
nmsFolS = dir([folderS '*.mat'])
matNmFlS = {nmsFolS.name}
matNmFlS = matNmFlS(1:end)

%% loop through all the files and put the arrays into a matrix ROCK
folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\rock\';
nmsFolR = dir([folderR '*.mat'])
matNmFlR = {nmsFolR.name}
matNmFlR = matNmFlR(1:end)
return

%% number of multimineralic grains with a given mineral (2 plot)

% for M = 2:length(mins)
    limMinAb = 0;
    close all
    f1 = figure
    M = 2
    minAb = mins(M)
    elMin = find(strcmp(mins,minAb))
    bins = 8

    subplot(1,2,2)
    for k = 1:length(matNmFlS)
        load([folderS matNmFlS{k}]);
        nRow = size(MnrlMtx,1);
        nCol = size(MnrlMtx,2);
        binC_mtx = [];
        mT_Mtx = [];
        count = 1;
        H = [];
        sz = [];

        [szM] = monoMin_by_mineral(MnrlMtx,elMin,limMinAb,nRow); %% here you call on a function to manipulate the MnrlMtx in whatever way you want
        [binC,yM] = bin_szHist(szM,bins);             
        [szA] = grainSizes_MnrlMtx(MnrlMtx,nRow);
        [binC,yA] = bin_szHist(szA,bins); 
        
        keyboard
        if mxS(k) == 1
            mfc = 'r';
        else
            mfc = [0.7 0.7 0.7];
        end
        if stS(k) == 1
            symt = '-';
        else
            symt = '--';
        end
    %     figure
    %     plot(sz,H,'.')
        hold on
        p(k)=plot((binC_mtx),mT_Mtx,symt,'linewidth',2,'color',mfc);   
        text(binC_mtx(8),mT_Mtx(8),labelS(k),'fontsize',12)

    end
    ylabel('Mean homogeneity')
    xlabel('log [Grain size (pixels)]')
    grid on
    legend([p(1) p(3) p(20) p(2)],{'MXS', 'MSS','MXNS','MSNS'},'location','southwest')
    % ylim([0.9 1])
    title(['Sediment - ' minAb{1}])
    % xlim([0 6])
    set(gca,'fontsize',16)
    sp2Ax = get(gca,'ylim')


    subplot(1,2,1)
    for k = 1:length(matNmFlR)
        load([folderR matNmFlR{k}]);
        nRow = size(MnrlMtx,1);
        nCol = size(MnrlMtx,2);
        binC_mtx = [];
        mT_Mtx = [];
        count = 1;
        H = [];
        sz = [];

        [H, sz] = homog_by_mineral(MnrlMtx,elMin,limMinAb,nRow);
        [binC_mtx,mT_Mtx] = bin_szA(sz,H,bins);
        if prR(k) == 1
            mfc = 'm';
        else
            mfc = [0.7 0.7 0.7]-0.4;
        end
        if stR(k) == 1
            symt = '-';
        else
            symt = '--';
        end
    %     figure
    %     plot(sz,H,'.')
        hold on
        p(k)=plot((binC_mtx),mT_Mtx,symt,'linewidth',2,'color',mfc);    
        text(binC_mtx(8),mT_Mtx(8),labelR(k),'fontsize',12)

    end
    ylabel('Mean homogeneity')
    xlabel('log [Grain size (pixels)]')
    grid on
    legend([p(1) p(4) p(18) p(2)],{'PS', 'MSS','PNS','MSNS'},'location','Southwest')
    % ylim([0.6 1])
    title(['Rock - ' minAb{1}])
    % xlim([0 6])
    set(gca,'fontsize',16)
    sp1Ax = get(gca,'ylim')
    
    ylims = [min([sp1Ax(1) sp2Ax(1)]), min([sp1Ax(2) sp2Ax(2)])];
    ylim(ylims);

    subplot(1,2,2)
    ylim(ylims)
%     savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\multiMineralByMineral\meanHvsLogGS_all_' minAb{1}])
    keyboard
% end

%% number of multimineralic grains with a given mineral (2 plot) All rock in one subplot, all seds in another

% limMinAbA_Arr = [45 95];
limMinAbA_Arr = [0];

% for M = 1:length(limMinAbA_Arr)
%     for M = 7
typeStr = 'Type_5c'
y_typeStr = 'Type 5c'
for L = 1:2
limStr = limMinAbA_Arr(L);
limStrS = num2str(limStr);
limMinAb = limStr/100;
txt = '$\rm \frac{Area\,of\,M\,in\,B}{Area\,of\,M\,in\,sample}$'
txt = 'mean homogeneity'
for M = 2:length(mins)
%     limMinAb = limMinAbA_Arr(M);
    
    close all
    f1 = figure('units','normalized','outerposition',[0 0 1 1])
%     M = 2
    minAb = mins(M)
    titleStr = [minAb{1} ' >' limStrS]
    elMin = find(strcmp(mins,minAb))
    bins = 8

    subplot(1,2,2)
    pK = 1
    mK = 1
    for k = 1:length(matNmFlS)
        load([folderS matNmFlS{k}]);
        nRow = size(MnrlMtx,1);
        nCol = size(MnrlMtx,2);
        binC_mtx = [];
        mT_Mtx = [];
        count = 1;
        H = [];
        sz = [];

        [H, szH] = homog_by_mineral(MnrlMtx,elMin,0,nRow);
        [binC,yH] = bin_szMean(szH,H,bins);
        
        if mxS(k) == 1
            mfc = 'r';
        else
            mfc = [0.7 0.7 0.7];
        end
        if stS(k) == 1
            symt = '-';
        else
            symt = '--';
        end
    %     figure
    %     plot(sz,H,'.')
        hold on
        y = yH;
%         y = 100*(yH./sum(yM));

        p(k)=plot(binC,y,symt,'linewidth',2,'color',mfc);   
        text(binC(bins),y(bins),labelS(k),'fontsize',12)
        yMtxS(k,:)=y;
    end
    yMeanSed(M,:) = nanmean(yMtxS);
    yMtxS_mxM(M,:) = nanmean(yMtxS(mxS,:));
    yMtxS_mxsM(M,:) = nanmean(yMtxS(logical(mxS.*stS),:));
    yMtxS_mxnsM(M,:) = nanmean(yMtxS(logical(mxS.*nstS),:));
    yMtxS_msM(M,:) = nanmean(yMtxS(msS,:));
    yMtxS_mssM(M,:) = nanmean(yMtxS(logical(msS.*stS),:));
    yMtxS_msnsM(M,:) = nanmean(yMtxS(logical(msS.*nstS),:));
    ylabel(['Percent - ' y_typeStr])
    xlabel('log [Grain size (pixels)]')
    grid on
    legend([p(1) p(3) p(20) p(2)],{'MXS', 'MSS','MXNS','MSNS'},'location','northwest')
    % ylim([0.9 1])
    title(['Sediment - ' titleStr])
    % xlim([0 6])
    set(gca,'fontsize',16)
    sp2Ax = get(gca,'ylim')
%     ylim([0 0.9])
    text(0.8,0.9,txt,'units','normalized','interpreter','latex','fontsize',16)

    subplot(1,2,1)
    pK = 1
    mK = 1
    for k = 1:length(matNmFlR)
%     for k = 3
        load([folderR matNmFlR{k}]);
        nRow = size(MnrlMtx,1);
        nCol = size(MnrlMtx,2);
        binC_mtx = [];
        mT_Mtx = [];
        count = 1;
        H = [];
        sz = [];
% 

        [H, szH] = homog_by_mineral(MnrlMtx,elMin,0,nRow);
        [binC,yH] = bin_szMean(szH,H,bins);
        
        if prR(k) == 1
            mfc = 'm';
        else
            mfc = [0.7 0.7 0.7]-0.4;
        end
        if stR(k) == 1
            symt = '-';
        else
            symt = '--';
        end
    %     figure
    %     plot(sz,H,'.')
        hold on
        y =yH;
%         y = 100*(yH./sum(yM));
        
        p(k)=plot(binC,y,symt,'linewidth',2,'color',mfc);   
        text(binC(bins),y(bins),labelR(k),'fontsize',12)
        yMtxR(k,:)=y;
    end
    yMeanRock(M,:) = nanmean(yMtxR);
    yMtxR_pM(M,:) = nanmean(yMtxR(prR,:));
    yMtxR_psM(M,:) = nanmean(yMtxR(logical(prR.*stR),:));
    yMtxR_pnsM(M,:) = nanmean(yMtxR(logical(prR.*nstR),:));
    yMtxR_msM(M,:) = nanmean(yMtxR(msR,:));
    yMtxR_mssM(M,:) = nanmean(yMtxR(logical(msR.*stR),:));
    yMtxR_msnsM(M,:) = nanmean(yMtxR(logical(msR.*nstR),:));
    ylabel(['Percent - ' y_typeStr])
    xlabel('log [Grain size (pixels)]')
    grid on
    legend([p(1) p(4) p(18) p(2)],{'PS', 'MSS','PNS','MSNS'},'location','northwest')
    % ylim([0.6 1])
    title(['Rock - ' titleStr])
    % xlim([0 6])
    set(gca,'fontsize',16)
    sp1Ax = get(gca,'ylim')
    
    ylims = [min([sp1Ax(1) sp2Ax(1)]), max([sp1Ax(2) sp2Ax(2)])];
    ylim(ylims);
%     ylim([0 0.9])

    subplot(1,2,2)
    ylim(ylims)
%     ylim([0 0.9])
    
    folderOut = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\' typeStr '\subplot_2\f' limStrS '\']
    mkdir(folderOut)
    savePDFfunction(f1,[folderOut minAb{1} '_1_' limStrS])
    saveJPGfunction(f1,[folderOut minAb{1} '_1_' limStrS])

% number of multimineralic grains with a given mineral (24 plot) (depends on running the 2 plot shown above)

    ylims(1)=0.6
    ylims(2)=0.8

    f2 = figure(2)

    for k = 1:length(SPIS)
        subplot(6,4,SPIS(k))

        if mxS(k) == 1
            mfc = 'r';
        else
            mfc = [0.7 0.7 0.7];
        end
        if SPIS(k) == 1 || SPIS(k) == 5 || SPIS(k) == 9 || SPIS(k) == 13 || SPIS(k) == 18 || SPIS(k) == 22
            ylabel('Percent')
        end
        if SPIS(k) == 13 || SPIS(k) == 22 || SPIS(k) == 19 || SPIS(k) == 20
            xlabel('log[ Grain size (pixels)]')
        end
        if SPIS(k) == 1 
            title(['MSNS - ' minAb{1} ' >' limStrS])
        end
        if SPIS(k) == 2 
            title('MSS')
        end
        if SPIS(k) == 3 
            title('MXNS')
        end
        if SPIS(k) == 4 
            title('MXS')
        end
        
        y = yMtxS(k,:)
        hold on
        pS(k) = plot(binC,y,'--','linewidth',2,'color',mfc);   
        text(binC(end),y(end),labelS(k),'fontsize',12)
        grid on
%         yl = get(gca,'ylim');
%         if yl(1)<ylims(1)
%             ylims(1)=yl(1);
%         end
%         if yl(2)>ylims(2)
%             ylims(2)=yl(2);
%         end
        
    end


    for k = 1:length(SPIR)
        subplot(6,4,SPIR(k))
       
        if prR(k) == 1
            mfc = 'm';
        else
            mfc = [0.7 0.7 0.7]-0.4;
        end
    %     figure
    %     plot(sz,H,'.')
        y = yMtxR(k,:)
        hold on
        pR(k) = plot(binC,y,'-','linewidth',2,'color',mfc);    
        text(binC(end),y(end),labelR(k),'fontsize',12)
%         yl = get(gca,'ylim');
%         if yl(1)<ylims(1)
%             ylims(1)=yl(1);
%         end
%         if yl(2)>ylims(2)
%             ylims(2)=yl(2);
%         end
    end

    for k = 1:length(SPIR)
        subplot(6,4,SPIR(k));
%         ylim([ylims]);
    end
    ax = subplot(6,4,23)
	ax.Visible = 'off'
    legend(ax,[pS(1), pS(2), pR(1), pR(2)],{'mixed seds','metased seds','plutonic rock','metased rock'},'location','southeast')
    
    ax = subplot(6,4,24)
	ax.Visible = 'off'
    text(0.1,0.7,y_typeStr,'fontsize',16)
    text(0.1,0,txt,'units','normalized','interpreter','latex','fontsize',16)
    
    folderOut = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\' typeStr '\subplot_24\f' limStrS '\']
    mkdir(folderOut)
    savePDFfunction(f2,[folderOut minAb{1} '_3_' limStrS])
    saveJPGfunction(f2,[folderOut minAb{1} '_3_' limStrS])

% in this section make a figure that has just the means from all types
    close all
    f3 = figure(3)
    subplot(2,2,1)
    hold on
    p1 = plot(binC,yMeanRock(M,:),'k-','linewidth',4)
    p2 = plot(binC,yMeanSed(M,:),'--','color','k','linewidth',4)
    p3 = plot(binC,yMtxR_pM(M,:),'-','color','m','linewidth',2)
    p6 = plot(binC,yMtxR_msM(M,:),'-','color',[0.3 0.3 0.3],'linewidth',2)
    p9 = plot(binC,yMtxS_mxM(M,:),'--','color','r','linewidth',2)
    p12 = plot(binC,yMtxS_msM(M,:),'--','color',[0.7 0.7 0.7],'linewidth',2)

    grid on
    legend([p1 p2 p3 p6 p9 p12],{'All rock','All seds',...
        'Plutonic rock','Metased. Rock',...
        'Mixed seds.','Metased. seds'},'fontsize',10,'location','best')
%     title([minAb{1} ' >' limStrS])
%     xlabel('log[ Grain size (pixels)]')
    ylabel(['Percent'])
    set(gca,'fontsize',16)
    
    subplot(2,2,3)
    hold on

    p7 = plot(binC,yMtxR_mssM(M,:),'-+','color',[0.3 0.3 0.3],'linewidth',2)
    p8 = plot(binC,yMtxR_msnsM(M,:),'-o','color',[0.3 0.3 0.3],'linewidth',2)

    grid on
    legend([p7 p8],{'Metased surge rock','Metased non-surge rock'},'fontsize',10,'location','best')
%     title([minAb{1} ' >' limStrS])
    xlabel('log[ Grain size (pixels)]')
    ylabel(['Percent'])
    set(gca,'fontsize',16)
    
    subplot(2,2,2)
    hold on

    p10 = plot(binC,yMtxS_mxsM(M,:),'--+','color','r','linewidth',2)
    p11 = plot(binC,yMtxS_mxnsM(M,:),'--o','color','r','linewidth',2)
    p13 = plot(binC,yMtxS_mssM(M,:),'--+','color',[0.7 0.7 0.7],'linewidth',2)
    p14 = plot(binC,yMtxS_msnsM(M,:),'--o','color',[0.7 0.7 0.7],'linewidth',2)
    grid on
    legend([p10 p11 p13 p14],{'MXS seds','MXNS seds','MSS seds','MSNS seds'},'fontsize',10,'location','best')
%     title([minAb{1} ' >' limStrS])
    xlabel('log[ Grain size (pixels)]')
%     ylabel(['Percent - ' y_typeStr])
    set(gca,'fontsize',16)
    
    ax = subplot(2,2,4)
	ax.Visible = 'off'
    text(0.1,0.1,[y_typeStr ', ' minAb{1} ' >' limStrS],'fontsize',16)
    text(0.1,0.3,txt,'units','normalized','interpreter','latex','fontsize',24)
    
    folderOut = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\' typeStr '\subplot_4\f' limStrS '\']
    mkdir(folderOut)
    savePDFfunction(f3,[folderOut minAb{1} '_2_' limStrS])
    saveJPGfunction(f3,[folderOut minAb{1} '_2_' limStrS])    
%     return
end
end

%%
close all
f1 = figure

% for i=1:length(limMinAbA_Arr)
for i=[2,6,9]
    
    subplot(1,2,1)
    hold on
%     plot(binC,yMeanRock(i,:),'k','linewidth',1)
%     text(binC(end),yMeanRock(i,end),num2str(limMinAbA_Arr(i)))
    plot(binC,yMtxR_pM(i,:),'m','linewidth',1)
    text(binC(end),yMtxR_pM(i,end),num2str(limMinAbA_Arr(i)))
    plot(binC,yMtxR_sM(i,:),'linewidth',1,'color',[0.3 0.3 0.3])
    text(binC(end),yMtxR_sM(i,end),num2str(limMinAbA_Arr(i)))
    
    subplot(1,2,2)
    hold on
%     plot(binC,yMeanSed(i,:),'linewidth',1,'color','k')
%     text(binC(end),yMeanSed(i,end),num2str(limMinAbA_Arr(i)))
    plot(binC,yMtxS_pM(i,:),'m','linewidth',1)
    text(binC(end),yMtxS_pM(i,end),num2str(limMinAbA_Arr(i)))
    plot(binC,yMtxS_sM(i,:),'linewidth',1,'color',[0.3 0.3 0.3])
    text(binC(end),yMtxS_sM(i,end),num2str(limMinAbA_Arr(i)))
end

subplot(1,2,1) 
grid on
ylabel('mean % for area % quartz > X')
xlabel('log[Grain size (pixels)]')
title('Rock')
ylim([0.05 0.55])
set(gca,'fontsize',16)

subplot(1,2,2) 
grid on
% ylabel('mean % for area % quartz > X')
xlabel('log[Grain size (pixels)]')
title('Sediment')
ylim([0.05 0.55])
set(gca,'fontsize',16)
% savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\percentGrainsOfMins\pGrainsvsLogGS_mean_BCT_' minAb{1}])



return
%% multimenarolgy of all grains all minerals combined (2plot)
close all
f1 = figure(1)
sp1 = subplot(1,2,1)
hold on
p=[]

stype = matNmFlR
bins = 8
mT_Mtx = zeros(length(stype),bins)        
for i = 1:length(stype)
%     subplot(6,4,SPIS(i))
    load([folderR stype{i}]);
    maxS(i) = max(sz_Arr);
    [binC_mtx(i,:),mT_Mtx(i,:)] = bin_szA(sz_Arr,T_Arr,bins);
    if prR(i) == 1
        mfc = 'm';
    else
        mfc = [0.7 0.7 0.7]-0.4;
    end
    if stR(i) == 1
        symt = '-';
    else
        symt = '--'
    end
    
    p(i) = plot(binC_mtx(i,:)*(25/34),mT_Mtx(i,:),symt,'linewidth',2,'color',mfc)
    text(binC_mtx(i,8)*(25/34),mT_Mtx(i,8),labelR(i),'fontsize',12)
end

ylabel('Mean homogeneity')
xlabel('log [Grain size (\mum)]')
grid on
legend([p(1) p(4) p(18) p(2)],{'PS', 'MSS','PNS','MSNS'},'location','Northwest')
ylim([0.6 1])
title('Rock')
xlim([0 6])
set(gca,'fontsize',16)

sp1 = subplot(1,2,2)
hold on

stype = matNmFlS
bins = 8
mT_Mtx = zeros(length(stype),bins)        
for i = 1:length(stype)-1
%     subplot(6,4,SPIS(i))
    load([folderS stype{i}]);
    maxS(i) = max(sz_Arr);
    [binC_mtx(i,:),mT_Mtx(i,:)] = bin_szA(sz_Arr,T_Arr,bins);
    if mxS(i) == 1
        mfc = 'r';
    else
        mfc = [0.7 0.7 0.7];
    end
    if stS(i) == 1
        symt = '-';
    else
        symt = '--'
    end
    hold on
    p(i) = plot(binC_mtx(i,:)*(25/34),mT_Mtx(i,:),symt,'linewidth',2,'color',mfc)
    text(binC_mtx(i,8)*(25/34),mT_Mtx(i,8),labelS(i),'fontsize',12)
end

ylabel('Mean homogeneity')
xlabel('log [Grain size (\mum)]')
grid on
legend([p(1) p(3) p(20) p(2)],{'MXS', 'MSS','MXNS','MSNS'},'location','northwest')
ylim([0.6 1])
title('Sediment')
xlim([0 6])
set(gca,'fontsize',16)

savePDFfunction(f1,'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\meanHvsLogGS_all')
%% multimineralogy of all minerals all samples (24 plot)
close all
f1 = figure(1)
hold on

bins = 8
mT_Mtx = zeros(length(matNmFlS),bins)        
for i = 1:length(SPIS)
    subplot(6,4,SPIS(i))
    load([folderS matNmFlS{i}]);
    maxS(i) = max(sz_Arr);
    [binC_mtx(i,:),mT_Mtx(i,:)] = bin_szA(sz_Arr,T_Arr,bins);
    if mxS(i) == 1
        mfc = 'r';
    else
        mfc = [0.7 0.7 0.7];
    end
    hold on
    plot(binC_mtx(i,:),mT_Mtx(i,:),'--','linewidth',2,'color',mfc)
    text(binC_mtx(i,end-1),mT_Mtx(i,end-1),labelS(i))
    if SPIS(i) == 1 || SPIS(i) == 5 || SPIS(i) == 9 || SPIS(i) == 13 || SPIS(i) == 18 || SPIS(i) == 22
        ylabel('H_m')
    end
    if SPIS(i) == 13 || SPIS(i) == 22 || SPIS(i) == 19 || SPIS(i) == 20
        xlabel('log[ Grain size (\mum)]')
    end
    ylim([0.6 0.95])
    xlim([0 6])
    grid on
    set(gca,'fontsize',12)

end

mT_Mtx = zeros(length(matNmFlR),bins)   
binC_mtx = zeros(length(matNmFlR),bins) 
for i = 1:length(SPIR)
    subplot(6,4,SPIR(i))
    load([folderR matNmFlR{i}]);
    maxS(i) = max(sz_Arr);
    [binC_mtx(i,:),mT_Mtx(i,:)] = bin_szA(sz_Arr,T_Arr,bins);
    if prR(i) == 1
        mfc = 'm';
    else
        mfc = [0.7 0.7 0.7]-0.4;
    end
    
    hold on
    plot(binC_mtx(i,:),mT_Mtx(i,:),'-','linewidth',2,'color',mfc)
    text(binC_mtx(i,end-1),mT_Mtx(i,end-1),labelR(i))

    ylim([0.6 0.95])
    xlim([0 6])
    grid on
    set(gca,'fontsize',12)

end

savePDFfunction(f1,'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\meanHvsLogGS_allSP_24')

return
%% histogram of samples for glacier 1
clear variables
close all

bins = 	8;
f1 = figure(1)

load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\seds\GL 01.mat')
[binC,yBin]=bin_szA(sz_Arr,T_Arr,bins)
hold on
% plot(sz_Arr,T_Arr,'k.')
p1 = plot((binC),yBin,'k-+','linewidth',2)
 
load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\rock\GL01_1.mat')
[binC,yBin]=bin_szA(sz_Arr,T_Arr,bins)
% plot(sz_Arr,T_Arr,'r.')
p2 = plot((binC),yBin,'r-+','linewidth',2)

load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\seds\GL1_1.mat')
[binC,yBin]=bin_szA(sz_Arr,T_Arr,bins)
% plot(sz_Arr,T_Arr,'b.')
p3 = plot((binC),yBin,'b-+','linewidth',2)

ylabel('Mean homogeneity')
xlabel('log [Grain size (pixels)]')
grid on
legend([p1 p3 p2],{'PG sed','BH sed', 'Rock'},'location','southeast')
% figure
% plot(lsa,T_Arr,'.')
% hold on
% plot(binC,yBin,'k-o','linewidth',2)

ylim([0.5 1])
title('Glacier 1')
set(gca,'fontsize',16)
% savePDFfunction(f1,'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\meanHvsLogGS_GL01')

return
%% this was for plotting glacier 1 samples

clear variables
BH = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\seds\GL1_1.mat')
RK = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\rock\GL01_1.mat')
SD = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\isolated_grains\seds\GL 01.mat')

f1 = figure

subplot(1,3,1)
p2 = plot(log(BH.sz_Arr),BH.T_Arr,'o','markerfacecolor',[0 0 1],'markeredgecolor','w','markersize',4)
grid minor
xlabel('log Grain size (pixels)')
ylabel('Homogeneity')
title('Borehole')
set(gca,'fontsize',16)

subplot(1,3,3)
p3 = plot(log(RK.sz_Arr),RK.T_Arr,'o','markerfacecolor',[1 0 0],'markeredgecolor','w','markersize',4)
xlabel('Grain size (pixels)')
grid minor
xlabel('log Grain size (pixels)')
title('Rock')
% title('Glacier 1 samples')
% set(gca,'fontsize',18)
set(gca,'fontsize',16)


subplot(1,3,2)
p1 = plot(log(SD.sz_Arr),SD.T_Arr,'o','markerfacecolor',[0 0 0],'markeredgecolor','w','markersize',4)
% legend([p1 p2 p3],{'PG Sed.','BH','Rock'},'location','southeast')
grid 
title('Sediment')
xlabel('log Grain size (pixels)')
set(gca,'fontsize',16)


savePDFfunction(f1,'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\HvsGS_GL1_3plot_log.pdf')

