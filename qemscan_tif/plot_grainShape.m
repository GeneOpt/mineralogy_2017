%%% in this script, I load the data generated from grainShape. I take the
%%% frequency spectra, which is composed of real and imaginary parts, and
%%% compute the mean spectra for groups based on the size and mineralogy.
%%% the minimum grian size spectra were computed on was 5^2 pixels

clear variables
run mineral_colors.m
%% loop through all the files and put the arrays into a matrix ROCK
folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\shape\rock\';
nmsFolR = dir([folderR '*.mat'])
matNmFlR = {nmsFolR.name}
matNmFlR = matNmFlR(1:end)

load([folderR matNmFlR{1}])

%%
close all
f1 = figure(1)
numP = size(mnrlMtxP,1);

elM = 1:numP;
gM = fSpec(elM,:);
realG = abs(gM);

for M = 2:length(mins)
    
    subplot(5,6,M-1)
    
    if M == 2
        for i = 1:size(realG,1)  
            y(i,:) = log(realG(i,:)/max(realG(i,:)));
            plot(y(i,:),'k.');
            hold on
        end
    end
    

    fracM = mnrlMtxP(:,M)./sum(mnrlMtxP,2);
    elM = find(fracM>0.4)
    ['M = ' num2str(M)]
    
    gM = fSpec(elM,:);
    realG = abs(gM);
    
    yM = []
    for i = 1:size(realG,1)
        hold on
        yM(i,:) = log(realG(i,:)/max(realG(i,:)));
        plot(yM(i,:),'.','color',min_col(M,:));
    end

    plot(mean(y),'k','linewidth',1)
    plot(mean(yM),'color','m','linewidth',2)
    grid on
%     ylabel('log(Normalized power)')
%     xlabel('Harmonic number (N)')
    ylim([-10 0])
    text(0.7,0.9,(minsN(M)),'units','normalized')

end
savePDFfunction(f1,'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\grain shape\plotTirty')

