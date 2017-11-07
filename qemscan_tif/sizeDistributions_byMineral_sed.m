clear all
close all
run loadSample_specs
run mineral_colors.m
    
xValR = [16 1 2 6 7 7 8 17 3 3 4 4 9 9 18 18 12 12 13 19 19 20 20 14 14 15]
xValS = [16 1 5 2 6 7 8 17 3 4 9 10 18 11 12 13 19 20 14 15 16]   
s2rA = [1 NaN; 2 NaN; NaN NaN; 3 NaN; 4 NaN; 5 6; 7 NaN; 8 NaN; 9 10; 11 12;...
    13 14; NaN NaN; 15 16; 16 NaN; 17 18; 19 NaN; 20 21; 22 23; 24 25; 26 NaN]
% get all the folder names (each folder is a rock sample)

fOut = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\fitGSD\islands\sed'

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
minP = 0.7420;
maxP = 62.0;
nBins = 10;

for i = 1:length(matNmS)

    mns = matNmS{i}
    varS = load([folderS mns]);
    
    mnrl_S = varS.mnrlMtx_c;
    [mnrlS,minsN,minNFull] = abbvMins(mnrl_S,5);
    islD = varS.isleD;
    
    fOutG = [fOut '\' matNmS{i}(1:5)];
    mkdir(fOutG)
    
    b = [];
    m = [];
    yLN_j_r = [];
    yLN_r = [];
    yW_j_r = [];
    yW_r = [];
    yFrac_r = [];
    sM = [];
    muM = [];
    lbdM = [];
    kM = [];
    sJ = [];
    muJ = [];
    lbdJ = [];
    kJ = [];
    mSL = [];
    bSL = [];

    for M = 1:length(fID)
%     for M = 13
        (fID{M})
        cnt = 1;
        elS = [];

%         for k = 1:size(mnrlS,1)
%             a = mnrlS(k,:);
%             aN = a(M)./sum(a(1:38));
%             if aN>0
%                 elS(cnt) = k;
%                 cnt = cnt+1;
%             end
% 
%         end   

%         gS = varS.D_c;
        islM = (islD.(fID{M}));
        x = islM;
%         x = gS;
        x(x==0) = [];

        [binC,binC_S,yBin,mBin,stdBin,sDat] = bin_szHist(x,x);
        
        y = yBin/sum(yBin);

        %% don't try to fit zeros or NaN's
        inY = isnan(y);
        elY = find(y==0);
        
        if sum(inY)==0
            if isempty(elY)==0
                elY = elY(1);
                y = y(1:elY-1);
                y = y/sum(y);
                binC = binC(1:elY-1);
            end
        end
        
        if length(y)>3 && sum(inY)==0
            
            %% try different fit types
            [yN,ySL,xSL,mSL(M),bSL(M)] = fitSL(x)
            [yW_j,kJ(M),lbdJ(M)] = fitW_J(y);
            [yLN_j,muJ(M),sJ(M)] = fitLN_J(y);

            [m(M), b(M), ~,~,~] = linReg(binC,log2(y));
            yLF = m(M)*binC+b(M);
            l2y = log2(y);

            binF = [binC,[binC(end)+0.5:0.5:-3.57]];
            binF = binC;

            parmhatW = wblfit(x);
            kM(M) = parmhatW(1);
            lbdM(M) = parmhatW(2);
            xW = linspace(minP,maxP,100);
            yW = wblpdf(2.^binF*1000,parmhatW(1),parmhatW(2));

            parmhatLN = lognfit(x);
            muM(M) = parmhatLN(1);
            sM(M) = parmhatLN(2);
            yLN = lognpdf(2.^binF*1000,parmhatLN(1),parmhatLN(2));

            %% do the plotting
            close all
            f1 = figure('Position',[0.0010    0.7610    1.9200    0.9673]*1e3);
            hold on
            p1 = plot(binC,y,'k-o','linewidth',3);
            yFrac = 2.^yLF;
            p2 = plot(binC,yFrac,'r-o','linewidth',1);
                yFrac_r(M) = sqrt(sum((y-yFrac).^2)/length(y));
            p4 = plot(binC,yW_j,'-o','linewidth',1,'color',[0 0.7 0]);
                yW_j_r(M) = sqrt(sum((y-yW_j).^2)/length(y));
            p6 = plot(binC,yLN_j,'-o','linewidth',1,'color',[0.9 0 0.9]);
                yLN_j_r(M) = sqrt(sum((y-yLN_j).^2)/length(y));
            p8 = plot(log2(xSL/1000),ySL,'b--s')
            p7 = plot(log2(xSL/1000),yN,'-s','linewidth',3,'color',[0.5 0.5 0.5])
                ySL_r(M) = sqrt(sum((yN-ySL).^2)/length(yN));

            grid on
            legend([p1 p2 p4 p6 p7 p8],{'Data (log2 binned)',['Fractal, rmd = ' num2str(yFrac_r(M),2)],...
                ['Weibull, rmd = ' num2str(yW_j_r(M),2)],...
                ['Log normal, rmd = ' num2str(yLN_j_r(M),2)],'Data (linear bin)',...
                ['Semi-log, rmd = ' num2str(ySL_r(M),2)]},...
                'location','Northeast','fontsize',16)
            title([matNmS{i}(1:5) ' - ' minsN{M} ' - Sediment'])
            ylabel('Number of grains (pdf)')
            xlabel('Grain size (-\phi)')
            set(gca,'fontsize',18)
    %         ylim([0 0.8])
            saveJPGfunction(f1,[fOutG '\' minsN{M}])
            savePDFfunction(f1,[fOutG '\' minsN{M}])

        else
            b(M) = 0;
            m(M) = 0;
            yLN_j_r(M) = NaN;
            yLN_r(M) = NaN;
            yW_j_r(M) = NaN;
            yW_r(M) = NaN;
            yFrac_r(M) = NaN;
            sM(M) = NaN;
            muM(M) = NaN;
            lbdM(M) = NaN;
            kM(M) = NaN;
            sJ(M) = NaN;
            muJ(M) = NaN;
            lbdJ(M) = NaN;
            kJ(M) = NaN;
            mSL(M)= NaN;
            bSL(M) = NaN;
            
        end

    end

paramM(:,:,i) = [m',b',yFrac_r',kJ',lbdJ',yW_j_r',muJ',sJ',yLN_j_r'];
 
end

save([fOut '\paramM_GSD_all.mat'],'paramM')
    










    
    