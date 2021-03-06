% this script is used to generate I and P plots for all glaciers in the
% islands\sizeDist folder. Each glacier has its own plot and is plotted
% with sediment and rock.

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


for i = 1:length(matNmS)

    mm = struct
    mnrlS = []
    mnrl_S = []
    gR = []
    mns = matNmS{i};
    varS = load([folderS mns])
    
    fn = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\islands\sizeDist\bins10\' matNmS{i}(1:5)]
    mkdir(fn)

    islS = varS.isleD;
    ptcS = varS.ptclD;
    gS = varS.D_c;
    mnrl_S = varS.mnrlMtx_c;
    [mnrlS,minsN,minNFull] = abbvMins(mnrl_S,5);
    
    r_i = s2rA(i,:);
    Hr = r_i(logical(1-isnan(r_i)));
    ct = 1
    if isempty(Hr)==0
        for j = 1:length(Hr)
            varR(j) = load([folderR matNmR{Hr(j)}])
            islR(j) = varR(j).isleD;
            ptcR(j) = varR(j).ptclD;
            [mm(j).R,minsN,minNFull] = abbvMins(varR(j).mnrlMtx_c,5);
            gR(j).R = varR(j).D_c;
        end
    end
   
    
%     for M = 1:length(fID)
    for M = 14
        close all
        p = [];
        
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
        
        f1 = figure
        p(1) = plot(binC,yBinI/sum(yBinI),'-o','color',colS(mxS(i)+1,:),'linewidth',2);
        hold on
        p(2) = plot(binC,yBinP/sum(yBinP),'--o','color',colS(mxS(i)+1,:),'linewidth',2);

        
        %% rock
        if isempty(Hr)==0
            for j = 1:length(Hr)
                
                mnrlR = mm(j).R;
                islR_M = islR(j).(fID{M});
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
              
                D_mR = gR(j).R(elR);
                [binC,binC_S,yBinI,mBin,stdBin,sDat] = bin_szHist(islR_M,islR_M);
                [binC,binC_S,yBinP,mBin,stdBin,sDat] = bin_szHist(D_mR,D_mR);
                if j == 1
                    cos = 0;
                elseif j == 2 & pr(Hr(j))+1 ~=2
                    cos = [0.3 0 0];
                elseif j == 2 & pr(Hr(j))+1 == 2
                    cos = 0;
                end
                    
                hold on
                p((2*j)+1) = plot(binC,yBinI/sum(yBinI),'-o','color',colR(pr(Hr(j))+1,:)+cos,'linewidth',2);
                hold on
                p((2*j)+2) = plot(binC,yBinP/sum(yBinP),'--^','color',colR(pr(Hr(j))+1,:)+cos,'linewidth',2);

                
            end           
        end
        
        
      
        
        grid on
        ylabel('Number of particles (pdf)');
        xlabel('Grain or island size (-\phi)');
        title([matNmS{i}(1:5) ' - ' minsN{M}]);
        if length(p)==2
            legend([p(1) p(2)],[matNmS{i}(1:5) ' (island)'],[matNmS{i}(1:5) ' (grain)'],...
                'location','northeast');
        end
        if length(p)==4
            legend([p(1) p(2) p(3) p(4)],[matNmS{i}(1:5) ' (island)'],[matNmS{i}(1:5) ' (grain)'],...
                [matNmR{r_i(1)}(1:8) '(island)'],[matNmR{r_i(1)}(1:8) '(grain)'],...
                'location','northeast');
        end
        if length(p)==6
            legend([p(1) p(2) p(3) p(4) p(5) p(6)],[matNmS{i}(1:5) ' (island)'],[matNmS{i}(1:5) ' (grain)'],...
                [matNmR{r_i(1)}(1:8) '(island)'],[matNmR{r_i(1)}(1:8) '(grain)'],...
                [matNmR{r_i(2)}(1:8) '(island)'],[matNmR{r_i(2)}(1:8) '(grain)'],...
                    'location','northeast');
        end
        ylim([0 0.8])
        set(gca,'fontsize',18);

        saveJPGfunction(f1,[fn '\' fID{M}])
        savePDFfunction(f1,[fn '\' fID{M}])
        [fn '\' fID{M}]


    end

    
end
