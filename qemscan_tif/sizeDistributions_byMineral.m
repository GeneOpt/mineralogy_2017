clear all
close all
run loadSample_specs
run mineral_colors.m
    
xValR = [16 1 2 6 7 7 8 17 3 3 4 4 9 9 18 18 12 12 13 19 19 20 20 14 14 15]
xValS = [16 1 5 2 6 7 8 17 3 4 9 10 18 11 12 13 19 20 14 15 16]   
s2rA = [1 NaN; 2 NaN; NaN NaN; 3 NaN; 4 NaN; 5 6; 7 NaN; 8 NaN; 9 10; 11 12;...
    13 14; NaN NaN; 15 16; 16 NaN; 17 18; 19 NaN; 20 21; 22 23; 24 25; 26 NaN]
% get all the folder names (each folder is a rock sample)

fOut = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\fitGSD\bins10M\sed'

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

    mns = matNmS{i};
    varS = load([folderS mns]);
    varB = load([folderS matNmS{i}]);
    
    mnrl_S = varS.mnrlMtx_c;
    [mnrlS,minsN,minNFull] = abbvMins(mnrl_S,5);
    islD = varS.isleD;
%     
%     fOutG = [fOut '\' matNmS{i}(1:5)];
    mkdir(fOut)
    
%     M = 1
%     for M = 1:length(fID)
% 
%         cnt = 1;
%         elS = [];
% 
%         for k = 1:size(mnrlS,1)
%             a = mnrlS(k,:);
%             aN = a(M)./sum(a(1:38));
%             if aN>0
%                 elS(cnt) = k;
%                 cnt = cnt+1;
%             end
% 
%         end   

        gS = varS.D_c;
%         islM = (islD.(fID{M}));
%         x = islM;
        x = gS;
        x(x==0) = [];
        
        [binC,binC_S,yBin,mBin,stdBin,sDat] = bin_szHist(x,x);
        
        y = yBin/sum(yBin);

        %% try different fit types

        [yW_j,kJ(i),lbdJ(i)] = fitW_J(y);
        [yLN_j,muJ(i),sJ(i)] = fitLN_J(y);
        
        [m(i), b(i), ~,~,~] = linReg(binC,log2(y));
        yLF = m(i)*binC+b(i);
        
        binF = [binC,[binC(end)+0.5:0.5:-3.57]];
        binF = binC;
        
        parmhatW = wblfit(x);
        kM(i) = parmhatW(1);
        lbdM(i) = parmhatW(2);
        xW = linspace(minP,maxP,100);
        yW = wblpdf(2.^binF*1000,parmhatW(1),parmhatW(2));
        
        parmhatLN = lognfit(x);
        muM(i) = parmhatLN(1);
        sM(i) = parmhatLN(2);
        yLN = lognpdf(2.^binF*1000,parmhatLN(1),parmhatLN(2));
        
        %% do the plotting
        close all
        f1 = figure('Position',[0.0010    0.7610    1.9200    0.9673]*1e3);
%         f1 = figure()
%         subplot(1,2,1)
        hold on
        p1 = plot(binC,y,'k-o','linewidth',3);
        yFrac = 2.^yLF;
        p2 = plot(binC,yFrac,'r-o','linewidth',2);
            yFrac_r(i) = sqrt(sum((y-yFrac).^2));
        p3 = plot(binF,yW,'-o','linewidth',2,'color',[0 0.4 0]);
            yW_r(i) = sqrt(sum((y-yW).^2));
        p4 = plot(binC,yW_j,'-o','linewidth',2,'color',[0 0.7 0]);
            yW_j_r(i) = sqrt(sum((y-yW_j).^2));
        p5 = plot(binF,yLN,'-o','linewidth',2,'color',[0.6 0 0.6]);
            yLN_r(i) = sqrt(sum((y-yLN).^2));
        p6 = plot(binC,yLN_j,'-o','linewidth',2,'color',[0.9 0 0.9]);
            yLN_j_r(i) = sqrt(sum((y-yLN_j).^2));
        grid on
        legend([p1 p2 p3 p4 p5 p6],{'Data',['Fractal, rmd = ' num2str(yFrac_r(i),2)],['Weibull (mle), rmd = ' num2str(yW_r(i),2)],...
            ['Weibull (rmm), rmd = ' num2str(yW_j_r(i),2)],['Log normal (mle), rmd = ' num2str(yLN_r(i),2)],...
            ['Log normal (rmm), rmd = ' num2str(yLN_j_r(i),2)]},...
            'location','Northeast','fontsize',16)
%         title([matNmS{i}(1:5) ' - ' minsN{M}])
        title([matNmS{i}(1:5) '- Sediment'])
        ylabel('Number of grains (pdf)')
        xlabel('Grain size (-\phi)')
        set(gca,'fontsize',18)
        ylim([0 0.8])
%         saveJPGfunction(f1,[fOut '\' matNmS{i}(1:5) '\' minsN{M}])
%         savePDFfunction(f1,[fOut '\' matNmS{i}(1:5) '\' minsN{M}])
        saveJPGfunction(f1,[fOut '\' matNmS{i}(1:5)])
        savePDFfunction(f1,[fOut '\' matNmS{i}(1:5)])
       
%     end  

end

paramM = [m',b',yFrac_r',kM',lbdM',yW_r',kJ',lbdJ',yW_j_r',muM',sM',yLN_r',muJ',sJ',yLN_j_r'];
save([fOut '\paramM_GSD_all.mat'],'paramM')

    










    
    