clear all
close all
run mineral_colors.m;
    
minP = 0.7424; 

cm = colormap('parula');
Hdiv = linspace(0,1,size(cm,1));

% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nms] = dir([folder '\*.mat']);
matNm = {nms.name}

% for H = 3:28  % for the rock
% for H = 1:length(matNm)    % for the sediment
for Hh = 1:length(matNm)   % for the sediment
   
    mn = matNm{Hh}
    nfN = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\islands\sed\' mn(1:5) '\']
    
    mkdir(nfN)
    
    for M = 1:45

        close all
        fname = dir([folder matNm{Hh}])
        load([folder matNm{Hh}])
        fieldN = fieldnames(isleD)

        y = isleD.(fieldN{M});
        x = ptclD.(fieldN{M});
        h = ptclH.(fieldN{M});

        thetaA = linspace(0,45.1,7);
        dTh = thetaA(2)-thetaA(1);
        binC = thetaA(1:end-1)+dTh/2
        for t = 1:length(thetaA)-1
            ll = thetaA(t);
            ul = thetaA(t+1);
            b_l = minP*(1-tand(ll));
            b_u = minP*(1-tand(ul));
            ly = x.*tand(ll)+b_l;
            uy = x.*tand(ul)+b_u;
            elB = find(y<uy & y>=ly);
            bin_pD = x(elB);
            bin_iD = y(elB);
            bin_h = h(elB);
            [bC,bC_S,yBin_P(t,:),mBin,stdBin,sDat] = bin_szHist(bin_pD,bin_pD);
            [bC,bC_S,yBin_I(t,:),mBin,stdBin,sDat] = bin_szHist(bin_iD,bin_iD);
            [bC,bC_S,yBin_pH(t,:),mBin_ph(t,:),stdBin,sDat] = bin_szHist(bin_pD,bin_h);
            [bC,bC_S,yBin_iH(t,:),mBin_ih(t,:),stdBin,sDat] = bin_szHist(bin_iD,bin_h);
            numB(t) = length(elB);
        end
        

        for i = 1:size(mBin_ph,1)
            for j = 1:size(mBin_ph,2)
                [hc,hcind_p(i,j)] = min(abs(Hdiv-mBin_ph(i,j)));
            end
        end

        for i = 1:size(mBin_ih,1)
            for j = 1:size(mBin_ih,2)
                [hc,hcind_i(i,j)] = min(abs(Hdiv-mBin_ih(i,j)));
            end
        end

%     %%
%     f1 = figure
%     yy = yBin_P/sum(sum(yBin_P))
%     b1 = bar3(yy)
%     
%     [nBar, nGroup] = size(yy);
%     
%     for i = 1:nBar
%         for j = 1:nGroup
%             text(j,i,yy(i,j)+0.01,num2str(mBin_ph(i,j),2))
%         end
%     end
%     
%     nColors  = size(get(gcf, 'colormap'), 1);
%     colorInd = hcind_p;
%     for i = 1:nGroup
%        c     = get(b1(i), 'CData');
%        color = repelem(repmat(colorInd(:, i), 1, 4), 6, 1);
%        set(b1(i), 'CData', color);
%     end
%     
%     set(gca,'XTickLabel',bC)
%     for i = 1:length(binC)
%         bcs{i} = num2str(binC(i),3)
%     end
%     set(gca,'YTickLabel',bcs)
%     set(gca,'view',[-225.7667   28.4000])
%     xlabel('Grain size (\phi)')
%     ylabel('Bin angle (^o)')
%     zlabel(['pdf (isl size) - ' fieldN{M} ', N = ' num2str(sum(sum(yBin_P)))])
%     set(gca,'fontsize',18)
%     colorbar
% %     savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\islands\GL 01\yBin_P_' mins_abbv{MAbb_i(M)}])
% %     saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\islands\GL 01\yBin_P_' mins_abbv{MAbb_i(M)}])
    %%
    f2 = figure
    yy = yBin_I/sum(sum(yBin_I))
    b1 = bar3(yy)
    
    [nBar, nGroup] = size(yy);
    
    for i = 1:nBar
        for j = 1:nGroup
            text(j,i,yy(i,j)+0.01,num2str(mBin_ih(i,j),2))
        end
    end
    
    nColors  = size(get(gcf, 'colormap'), 1);
    colorInd = hcind_i;
    for i = 1:nGroup
       c     = get(b1(i), 'CData');
       color = repelem(repmat(colorInd(:, i), 1, 4), 6, 1);
       set(b1(i), 'CData', color);
    end
    
    set(gca,'XTickLabel',bC)
    for i = 1:length(binC)
        bcs{i} = num2str(binC(i),3);
    end
    set(gca,'YTickLabel',bcs)
    set(gca,'view',[-225.7667   28.4000])
    xlabel('Grain size (\phi)')
    ylabel('Bin angle (^o)')
    zlabel(['pdf (isl size) - ' fieldN{M} ', N = ' num2str(sum(sum(yBin_I)))])
    set(gca,'fontsize',18)
    savePDFfunction(f2,[nfN 'yBin_I_' fieldN{M}])
    saveJPGfunction(f2,[nfN 'yBin_I_' fieldN{M}])
    %%
    f3 = figure
    bar(binC,numB/sum(numB))
    ylabel('pdf')
    xlabel('Bin angle (^o)')
    grid on
    title(fieldN{M})
    set(gca,'fontsize',18)

    savePDFfunction(f3,[nfN 'numB_' fieldN{M}])
    saveJPGfunction(f3,[nfN 'numB_' fieldN{M}])
%%
    f4 = figure
    for i = 1:length(h)
        [hc,hcind_a(i)] = min(abs(Hdiv-h(i)));
        hold on
        plot(x(i),y(i),'.','color',cm(hcind_a(i),:))
    end
    grid on
    ylabel('Grain diameter island (\mum)')
    xlabel('Grain diameter particle (\mum)')
    rl=refline(1,0)
    rl.Color = 'k';
    title(fieldN{M})
    set(gca,'fontsize',18)
    
    savePDFfunction(f4,[nfN 'prtcl_isl_' fieldN{M}])
    saveJPGfunction(f4,[nfN 'prtcl_isl_' fieldN{M}])

    end
end
