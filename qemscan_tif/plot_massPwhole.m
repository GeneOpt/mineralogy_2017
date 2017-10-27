% clear variables
close all

load('massPwhole.mat')
run mineral_colors

% massProck(:,39) = massProck(:,4)+massProck(:,5)+massProck(:,6)+massProck(:,7)
% massProck(:,40) = massProck(:,9)+massProck(:,10)+massProck(:,11)+massProck(:,12)
% massProck(:,41) = massProck(:,14)+massProck(:,15)
% massProck(:,42) = massProck(:,18)+massProck(:,19)
% massProck(:,43) = massProck(:,21)+massProck(:,22)
% massProck(:,44) = massProck(:,27)+massProck(:,28)
% 
% massPsed(:,39) = massPsed(:,4)+massPsed(:,5)+massPsed(:,6)+massPsed(:,7)
% massPsed(:,40) = massPsed(:,9)+massPsed(:,10)+massPsed(:,11)+massPsed(:,12)
% massPsed(:,41) = massPsed(:,14)+massPsed(:,15)
% massPsed(:,42) = massPsed(:,18)+massPsed(:,19)
% massPsed(:,43) = massPsed(:,21)+massPsed(:,22)
% massPsed(:,44) = massPsed(:,27)+massPsed(:,28)
% 

% massPsed(:,1) = []
% massProck(:,1) = []
% 

% save('massPwhole.mat','massPsed','massProck')


[datR,varsR] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',1,'B1:I27');
[datS,varsS] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',2,'B1:I27');

minsN = {'Qtz','Ksp','Ab','An25','AN50','An75','Mscv','Bti low','Bti Mg','Bti int',...
    'Bti Fe','Kao','Chl Fe','Chl Mg','Mg Clay','Mg Si','Ill Smec','Ill Smec Fe','Cal',...
    'Dol','Dol Fe','Fe ox si','Prt','Gyp Anh','Hlt','Rtl Ilm','Ilm','Ttn','Lmt','Cpx',...
    'Fe Amph','Epd Zo','Apt','Trml','Zrc','Calg','AlOx',...
    'Plag','Biot','Chlor','Ill smec c','Dol c','Ilmen'};

minNFull = {'Quartz', 'K Feldspar', 'Plagioclase Ab', 'Plagioclase An25',...
'Plagioclase An50', 'Plagioclase An75','Muscovite', 'Biotite (Low Fe & Mg)',...
'Biotite (Mg-rich)', 'Biotite (Intermediate)', 'Biotite (Fe-rich)', 'Kaolinite',...
'Fe Chlorite', 'Mg Chlorite', 'Mg Clays', 'Mg Silicate', 'Illite & illite-smectite',...
'Fe-Illite & illite-smectite','Calcite', 'Dolomite', 'Ferroan Dolomite',...
'Fe Oxide & siderite', 'Pyrite', 'Gypsum / Anhydrite', 'Halite', 'Rutile & Ilmenite',...
'Ilmenite', 'Titanite', 'Laumontite', 'Clinopyroxene', 'Fe Amphibole', 'Epidote / Zoisite',...
'Apatite', 'Tourmaline', 'Zircon','Calgon', 'Aluminium Oxide',...
'Plagioclase','Biotite','Chlorite','Illite Smectite c','Dolomite','Ilmenite'};
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
xValR = [16 1 2 6 7 7 8 17 3 3 4 4 9 9 18 18 12 12 13 19 19 20 20 14 14 15]

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
xValS = [16 1 5 2 6 7 8 17 3 4 9 10 18 11 12 13 19 20 14 15 16]   


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


%% plot individual ratios

elR = find(strcmp(samplesR,'GL01_1'))
elS = find(strcmp(samplesS,'GL 01'))
elM2 = find(strcmp(minsN,'Chl Fe'))
elM1 = find(strcmp(minsN,'Ab'))

rR = massProck(elR,elM1)/massProck(elR,elM2)
rS = massPsed(elS,elM1)/massProck(elS,elM2)

%% in this plot you can plot the glacier on x axis and the the rock and sed mass percent both 
% show up on the y axis
close all

% justClay = [7:15,17,18];
% massPs = massPsed(:,justClay);
% mpsN = massPs./repmat(sum(massPs,2),[1,length(justClay)]);
% massPr = massProck(:,justClay)
% mprN = massPr./repmat(sum(massPr,2),[1,length(justClay)]);

mpsN = massPsed(:,elM1)./massPsed(:,elM2);
mprN = massProck(:,elM1)./massProck(:,elM2);

% mnrl = 'Qtz'
% for j = 1:length(minsN)
for j = 1
    f1 = figure
    subplot(2,1,1)
 

    % elM = find(strcmp(mnrl,minsN))
    M_s = mpsN(:,j)


    colA = [0.7 0.7 0.7;1 0.7 0.7]
    mt = {'^','o'}

    for i = 1:length(xValS)

        hold on
        plot(xValS(i),M_s(i),['k' mt{mxS(i)+1}],'markerfacecolor', colA(stS(i)+1,:))
        text(xValS(i),M_s(i),labelS{i})

    end

    hold on
    M_r = mprN(:,j)

    colA = [0.1 0.1 0.1;1 0 0]
    mt = {'o','^'}

    for i = 1:length(xValR)

        hold on
        plot(xValR(i),M_r(i),['k' mt{msR(i)+1}],'markerfacecolor', colA(stR(i)+1,:))
        text(xValR(i),M_r(i),labelR{i})

    end
    xlim([0 22])
    grid on
    ylabel('Ratio')
    title([minsN{elM1} ' / ' minsN{elM2}])
    % XX = 1:length(labelS)
    set(gca(),'XTickLabel',[])
    set(gca,'fontsize',18)
    % rotateXLabels(gca(),60)


%% in this plot you can make a scatter plot of the difference between rock and sed. 

    % mnrl = 'Qtz'
    subplot(2,1,2)
%     for j = 1:length(minsN)


        % elM = find(strcmp(mnrl,minsN))
%         M_s = massPsed(:,j)
% 
%         M_r = massProck(:,j)
% 
        colA = [0.3 0.3 0.3;1 0.2 0.2]
%         mt = {'o','^'}
        diff_A = []
        for i = 1:length(xValR)

            xV = xValR(i)
            el = find(xValS==xV)
            el = el(1)
            M_si = M_s(el)
            M_diff = M_si-M_r(i)

            hold on
            plot(xValR(i),M_diff,['k' mt{msR(i)+1}],'markerfacecolor', colA(stR(i)+1,:))
            text(xValR(i),M_diff,labelR{i})
            M_diffA(i) = M_diff

        end
%         M_diffA(2) = []
        mDA = nanmean(M_diffA)
        [h,p] = ttest(M_diffA)
        text(0.1,0.1,['\mu = ' num2str(mDA,2) ', p = ' num2str(p,2)],'units','normalized')

        xlim([0 22])
        grid on 
%         title(minsN{j})
        ylabel('difference in ratio (sediment - rock)')
        set(gca,'fontsize',18)
        set(gca(),'XTickLabel',[])

        savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\massPercent\ratios\' minsN{elM1} '_' minsN{elM2}])
        saveJPGfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\massPercent\ratios\' minsN{elM1} '_' minsN{elM2}])
    %%
        % XX = 1:length(labelS)
        % set(gca(),'XTick',XX, 'XTickLabel',labelS)
        % rotateXLabels(gca(),60)
% close all
end




