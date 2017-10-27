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
    13 14; NaN NaN; 15 16; 15 NaN; 17 18; 19 NaN; 20 21; 22 23; 24 25; 26 NaN]
% get all the folder names (each folder is a rock sample)

folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\concatenated files\';
[nmsR] = dir([folderR '\*.mat']);
matNmR = {nmsR.name}

folderS = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\concatenated files\';
[nmsS] = dir([folderS '\*.mat']);
matNmS = {nmsS.name}
% for H = 3:28  % for the rock
% for H = 1:length(matNm)    % for the sediment
for Hs = 19


    mns = matNmS{Hs};
    varS = load([folderS mns])
    fNm = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\mineralAssociations\together\' mns(1:5) '\'];
    mkdir(fNm)
    
    MAM_cS = varS.MAM_c;
    MAM_all_cS = varS.MAM_all_c;

    AbM = 0

    if AbM == 1
        [MAM_c,minsN,minNFull] = abbvMins(MAM_c,4)
        [MAM_all_c] = abbvMins(MAM_all_c,3)
    end
    
%%
    for P = 1:size(MAM_all_cS,1)
        
        close all
        f1 = figure
        tArrS = MAM_all_cS(P,:);
        [tArrS,isr] = sort(tArrS,'descend');
        numPtS = sum(tArrS);
        yS = log10(tArrS/numPtS*100);
        h1 = plot(yS,'k-^') ;

        r_i = s2rA(Hs,:);
        Hr = r_i(logical(1-isnan(r_i)));
        ct = 1
        if isempty(Hr)==0
            for Hr_i = 1:length(Hr)
                      
                mnR = matNmR{Hr(Hr_i)};
                varR = load([folderR mnR]);
                MAM_cR = varR.MAM_c;
                MAM_all_cR = varR.MAM_all_c;

                tArrR = MAM_all_cR(P,:);
                tArrR= tArrR(isr);
                numPtR = sum(tArrR);
                yR = log10(tArrR/numPtR*100);
                hold on
            
                if pr(Hr(Hr_i)) == 1
                    text(3,yR(3),labelR(Hr(Hr_i)))
                    symb = '-^'
                    colr = 'm'

                elseif pr(Hr(Hr_i)) == 0
                    text(5,yR(5),labelR(Hr(Hr_i))) 
                    
                    colr = [0.3 0.3 0.3]
                    if ct == 1
                        symb = '-v'
                    else
                        symb = '--^'
                    end
                    ct = ct+1
                    ct = ct+1
                end
                h2 = plot(yR,symb,'color',colr);
                legend([h1, h2],{'Sediment','Rock'})
            end
        end
        
        grid on
        ylabel(minNFull{P})
        set(gca,'xtick',1:length(minsN),'xticklabel',minsN(isr))
        rotateXLabels(gca,60)
        set(gca,'fontsize',17) 
        
            
        saveJPGfunction(f1,[fNm mins{P}])
        savePDFfunction(f1,[fNm mins{P}])

    end
    

end