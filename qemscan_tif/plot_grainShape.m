%%% in this script, I load the data generated from grainShape. I take the
%%% frequency spectra, which is composed of real and imaginary parts, and
%%% compute the mean spectra for groups based on the size and mineralogy.
%%% the minimum grian size spectra were computed on was 5^2 pixels

clear variables
run mineral_colors.m
close all

%% loop through all the files and put the arrays into a matrix ROCK
folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\shape\rock\';
nmsFolR = dir([folderR '*.mat'])
matNmFlR = {nmsFolR.name}
matNmFlR = matNmFlR(1:end)

fSpecM = [];
mnrlMtxPM = [];

for H = 1:length(matNmFlR)
    
    load([folderR matNmFlR{H}])

    fSpecM = [fSpecM;fSpec];
    mnrlMtxPM = [mnrlMtxPM;mnrlMtxP];

end
%%

f1 = figure(1)
numP = size(mnrlMtxPM,1);

fSpecM = fSpecM(:,1:20);

elM = 1:numP;
gM = fSpecM(elM,:);
realG = abs(gM);

WHATMIN = 7

for M = WHATMIN
    
%     subplot(5,6,M-1)
    hold on
    
%     if M == 2
        for i = 1:size(realG,1)  
            y(i,:) = log(realG(i,:)/max(realG(i,:)));
%             plot(y(i,:),'k.');
            hold on
        end
%     end
    

    fracM = mnrlMtxPM(:,M)./sum(mnrlMtxPM,2);
    elM = find(fracM>0.4)
    ['M = ' num2str(M)]
    
    gM = fSpecM(elM,:);
    realG = abs(gM);
    
    yM = []
    for i = 1:size(realG,1)
        hold on
        yM(i,:) = log(realG(i,:)/max(realG(i,:)));
%         plot(yM(i,:),'o','color',min_col(M,:),'markerfacecolor',min_col(M,:),'markersize',1);
    end

    h1 = plot(mean(y),'k','linewidth',1)
    h2 = plot(mean(yM),'color','m','linewidth',2)
    grid on
%     ylabel('log(Normalized power)')
%     xlabel('Harmonic number (N)')
    ylim([-10 0])
%     text(0.7,0.9,(minsN(M)),'units','normalized')
  

end

%% loop through all the files and put the arrays into a matrix ROCK
folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\shape\sediment\';
nmsFolR = dir([folderR '*.mat'])
matNmFlR = {nmsFolR.name}
matNmFlR = matNmFlR(1:end)

fSpecM = [];
mnrlMtxPM = [];

for H = 1:length(matNmFlR)
    
    load([folderR matNmFlR{H}])

    fSpecM = [fSpecM;fSpec];
    mnrlMtxPM = [mnrlMtxPM;mnrlMtxP];

end
%%

f1 = figure(1)
numP = size(mnrlMtxPM,1);

fSpecM = fSpecM(:,1:20);

elM = 1:numP;
gM = fSpecM(elM,:);
realG = abs(gM);

for M = WHATMIN
    
%     subplot(5,6,M-1)
    hold on
    
%     if M == 2
        for i = 1:size(realG,1)  
            y(i,:) = log(realG(i,:)/max(realG(i,:)));
%             plot(y(i,:),'k.');
            hold on
        end
%     end
    

    fracM = mnrlMtxPM(:,M)./sum(mnrlMtxPM,2);
    elM = find(fracM>0.4)
    ['M = ' num2str(M)]
    
    gM = fSpecM(elM,:);
    realG = abs(gM);
    
    yM = []
    for i = 1:size(realG,1)
        hold on
        yM(i,:) = log(realG(i,:)/max(realG(i,:)));
%         plot(yM(i,:),'k^','markerfacecolor',min_col(M,:),'markersize',1);
    end

    h3 = plot(mean(y),'k--','linewidth',1)
    h4 = plot(mean(yM),'--','color','m','linewidth',2)
    grid on
%     ylabel('log(Normalized power)')
%     xlabel('Harmonic number (N)')
    ylim([-10 0])
    text(0.7,0.9,(minsN(M)),'units','normalized')

end

legend([h1 h2 h3 h4],{'Rock - mean of all minerals','Rock - mean of muscovite','Sediment - mean of all minerals','Sediment - mean of muscovite'},'location','northeast')
title('Glacier 02 muscovite')
xlabel('Harmonic number')
ylabel('log(normalized power)')
set(gca,'fontsize',18)
ylim([-6 -2])
xlim([1 20])
saveJPGfunction(f1,'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\grain shape\NWG2017\msc_rockSed')




