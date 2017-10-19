close all
clear all
run loadSample_specs

varS = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\GL 01.mat')
varR = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\GL01RS01.mat')
varB = load('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\GL01B.mat')

% for M = 2:length(mins)
for M = [4,6,12,14,19]

    D_s = varS.D_c;
    D_s_nPhi = log2(D_s/1000);
    D_r = varR.D_c;
    D_r_nPhi = log2(D_r/1000);
    D_b = varB.D_c;
    D_b_nPhi = log2(D_b/1000);
    m_s = varS.mnrlMtx_c;
    [m_s,minsN_abbv,minNFull_abbv] = abbvMins(m_s,2);
    H_s = varS.H;
    H_s_ab = varS.H_ab;
    H_r = varR.H;
    m_r = varR.mnrlMtx_c;
    [m_r,minsN_abbv,minNFull_abbv] = abbvMins(m_r,2);
    H_r_ab = varR.H_ab;
    H_b = varB.H;
    H_b_ab = varB.H_ab;
    m_b = varB.mnrlMtx_c;
    [m_b,minsN_abbv,minNFull_abbv] = abbvMins(m_b,2);

    ByMin = 1
    if ByMin == 1
        norm_M_s = m_s./repmat(sum(m_s,2),[1, size(m_s,2)]);
        elGt = find(norm_M_s(:,M)>0.05);
        H_s = H_s(elGt);
        D_s = D_s(elGt);

        norm_M_r = m_r./repmat(sum(m_r,2),[1, size(m_r,2)]);
        elGt = find(norm_M_r(:,M)>0.05);
        H_r = H_r(elGt);
        D_r = D_r(elGt);

        norm_M_b = m_b./repmat(sum(m_b,2),[1, size(m_b,2)]);
        elGt = find(norm_M_b(:,M)>0.05);
        H_b = H_b(elGt);
        D_b = D_b(elGt);    
    end
    %% plot the homogeneity for all minerals
    close all
    f1 = figure(10)
    [binC,binC_S,yBinS,mBinS,stdBinS,sDatS] = bin_szHist(D_s,H_s);
    % plot(D_s_nPhi,H_s,'k.')
    hold on
    p1 = plot(binC,mBinS,'ko','markerfacecolor','k','markersize',10);

    [binC,binC_R,yBinR,mBinR,stdBinR,sDatR] = bin_szHist(D_r,H_r);
    % plot(D_s_nPhi,H_s,'r.','markersize',3)
    hold on
    p3 = plot(binC,mBinR,'ro','markerfacecolor','r','markersize',10);


    [binC,binC_B,yBinB,mBinB,stdBinB,sDatB] = bin_szHist(D_b,H_b);
    % plot(D_s_nPhi,H_s,'b.','markersize',2)
    hold on
    p2 = plot(binC,mBinB,'bo','markerfacecolor','b','markersize',10);
    % f1 = fill([binC binC(end:-1:1)],[mBinS-stdBinS mBinS(end:-1:1)+stdBinS(end:-1:1)],'b')
    % f1.FaceAlpha = 0.2
    % f1.LineStyle = 'none'
    grid on
    ylabel('Mean homogeneity')
    xlabel('Grain size (-\phi)')
    legend([p1 p2 p3],{'PG Sediment','BH Sediment','Rock'},'location','southeast');
    set(gca,'fontsize',20)
    xlim([-11 -5]);


    %% here you can do the stats to see if they differ
    cAll = [];

    % between sed 
    for i = 1:length(binC)-1
        dat1 = sDatS.bin{i};
        dat2 = sDatS.bin{i+1};
        y = [dat1,dat2];
        bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
        p = kruskalwallis(y,bina,'off');
        cAll(i) = p;
    end

    elP = find(cAll<0.01);
    figure(10)
    for i = 1:length(elP)
        hold on
        plot([binC(elP(i)) binC(elP(i)+1)],[mBinS(elP(i)) mBinS(elP(i)+1)],'k','linewidth',2)
    end

    %% for the rock
    cAll = []

    for i = 1:length(binC)-1
        dat1 = sDatR.bin{i};
        dat2 = sDatR.bin{i+1};
        y = [dat1,dat2];
        bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
        p = kruskalwallis(y,bina,'off');
        cAll(i) = p;
    end

    elP = find(cAll<0.01);
    figure(10)
    for i = 1:length(elP)
        hold on
        plot([binC(elP(i)) binC(elP(i)+1)],[mBinR(elP(i)) mBinR(elP(i)+1)],'r','linewidth',2)
    end

    %% for the BH
    cAll = [];

    for i = 1:length(binC)-1
        dat1 = sDatB.bin{i};
        dat2 = sDatB.bin{i+1};
        y = [dat1,dat2];
        bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
        p = kruskalwallis(y,bina,'off');
        cAll(i) = p;
    end

    elP = find(cAll<0.01);
    figure(10)
    for i = 1:length(elP)
        hold on
        plot([binC(elP(i)) binC(elP(i)+1)],[mBinB(elP(i)) mBinB(elP(i)+1)],'b','linewidth',2);
    end


    %% for the sed vs rock
    cAll = [];

    for i = 1:length(binC)
        dat1 = sDatR.bin{i};
        dat2 = sDatS.bin{i};
        y = [dat1,dat2];
        bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
        p = kruskalwallis(y,bina,'off');
        cAll(i) = p;
    end

    elP = find(cAll<0.01);
    figure(10)
    for i = 1:length(elP)
        hold on
        plot([binC(elP(i)) binC(elP(i))],[mBinR(elP(i)) mBinS(elP(i))],'color',[0.5 0.5 0.5],'linewidth',2);
    end


    %% for the BH vs rock
    cAll = [];

    for i = 1:length(binC)
        dat1 = sDatR.bin{i};
        dat2 = sDatB.bin{i};
        y = [dat1,dat2];
        bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
        p = kruskalwallis(y,bina,'off')
        cAll(i) = p;
    end

    elP = find(cAll<0.01);
    figure(10)
    for i = 1:length(elP)
        hold on
        plot([binC(elP(i))-0.1 binC(elP(i))-0.1],[mBinR(elP(i)) mBinB(elP(i))],'color',[0.5 0.5 0.5],'linewidth',2);
    end

    %% for the sedsed vs BH
    cAll = [];

    for i = 1:length(binC)
        dat1 = sDatB.bin{i};
        dat2 = sDatS.bin{i};
        y = [dat1,dat2];
        bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
        p = kruskalwallis(y,bina,'off')
        cAll(i) = p;
    end

    figure(10)
    elP = find(cAll<0.01);
    for i = 1:length(elP)
        hold on
        plot([binC(elP(i))+0.1 binC(elP(i))+0.1],[mBinB(elP(i)) mBinS(elP(i))],'color',[0.5 0.5 0.5],'linewidth',2);
    end

    %% and finally with the mean of all minerals
    xlms = get(gca,'xlim')
    hold on
    plot(xlms,[mean(H_s) mean(H_s)],'k')
    hold on
    plot(xlms,[mean(H_r) mean(H_r)],'r')
    hold on
    plot(xlms,[mean(H_b) mean(H_b)],'b')
        
    dat1 = H_s;
    dat2 = H_r;
    y = [dat1,dat2];
    bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
    p = kruskalwallis(y,bina,'off')
    p01 = p

    dat1 = H_s;
    dat2 = H_b;
    y = [dat1,dat2];
    bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
    p = kruskalwallis(y,bina,'off')
    p02 = p
    
    dat1 = H_r;
    dat2 = H_b;
    y = [dat1,dat2];
    bina = [ones(1,length(dat1)),zeros(1,length(dat2))];
    p = kruskalwallis(y,bina,'off')
    p03 = p
    figure(10)
    text(0.05,0.15,['All sizes \newline H_{s-r}, p=' num2str(p01,2) '\newline H_{s-b}, p=' num2str(p02,2) '\newline H_{r-b}, p=' num2str(p03,2)],'fontsize',16,'units','normalized')
 %%   
%     edges = [0:.1:1]
%     binC = edges(1:end-1)+0.05
%     hS = histogram(H_s,edges)
%     hsv = hS.Values
%     hR = histogram(H_r,edges)
%     hsr = hR.Values
%     hB = histogram(H_b,edges)
%     hsb = hB.Values
%     yB = [hsv/sum(hsv);hsr/sum(hsr);hsb/sum(hsb)];
% 
%     figure
%     bar(binC,yB')
%     
    title([minNFull_abbv{M} ' - combined']);
    savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\homogeneity\GL 01\pdf\meanH_' minsN_abbv{M} '_c'])
    saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\homogeneity\GL 01\jpg\meanH_' minsN_abbv{M} '_c'])

end




