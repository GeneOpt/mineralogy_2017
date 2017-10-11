%% read in the qemscan data
% clear all
% run read_in_qemscan
% 
% save('tempDatDump.mat','var_ss','var_BMMP_F','dat_BMMP_F','varMD','datMD')

%%
clear variables
%7,9. rock was from a different location that mapping
P10 = [4.64 26.29 NaN 3.22 7.97 NaN 7.5 5.15 9.43 9.43 NaN NaN 8.22 8.22 ...
    NaN 5.8 19.06 NaN 10.37 5.49 NaN 5.46 NaN NaN 6.06 NaN]


load tempDatDump.mat
%% plotting specs

mt_pr = 'ko'
mt_ms = 'ko'
mfc_s = [1 0 0]
mfc_ns = [0.7 0.7 0.7]
ms = 2

%% samples 
% samples = {'12/25','GL02RS02','GL04RS03','GL05RS02','GL06RS02','GL06RS04',...
%     'GL07RS02','GL08RS01','GL09RS01A','GL09RS01B','GL10RS01','GL10RS02',...
%     'GL11RS01','GL11RS02','GL14RS03','GL14RS02','GL16RS1','GL16RS02',...
%     'GL17RS01','GL18RS01','GL18RS02','GL19RS01','GL19RS02','GL20RS1',...
%     'GL20RS02','GL21RS01'} %this is the full list
% 
samples = {'12/25','GL02RS02','GL04RS03','GL05RS02','GL06RS02','GL06RS04',...
    'GL07RS02','GL08RS01','GL09RS01A','GL09RS01B','GL10RS01','GL10RS02',...
    'GL11RS01','GL11RS02','GL14RS03','GL14RS02','GL16RS02',...
    'GL17RS01','GL18RS01','GL18RS02','GL19RS02','GL20RS1'} %%this is maybe a modified list

for i = 1:length(samples)
    keepLab(i) = find(strcmp(var_ss(:,1),samples(i)))
end

st = strcmp(var_ss(keepLab,3),'S')
nst = strcmp(var_ss(keepLab,3),'NS')
gt = var_ss(2:end,3)
p_r = strcmp(var_ss(keepLab,4),'P')
Gb_r = strcmp(var_ss(keepLab,4),'Gb')
ms_r = strcmp(var_ss(keepLab,4),'MS')
rt = var_ss(keepLab,4)
mx_c = strcmp(var_ss(keepLab,5),'MX')
ms_c = strcmp(var_ss(keepLab,5),'MS')
group = var_ss(keepLab,6)
lab = {'1','2','4','5','6a','6b','7','8','9a','9b','10a','10b',...
    '11a','11b','14a','14b','16a','16b','17','18a','18b','19a','19b',...
    '20a','20b','21'}
lab = lab(keepLab-1)

%% minerals
mins = {'Qtz','Ksp','An0','An25','An50','An75','Mscv','Bti_low','Bti_Mg','Bti_int',...
    'Bti_Fe','Kao','Chl_Fe','Chl_Mg','Mg_Clay','Mg_Si','Ill_Smec','Ill_Smec_Fe','Cal',...
    'Dol','Dol_Fe','Fe_Ox_Sid','Prt','Gyp_Anh','Hlt','Rtl_Ilm','Ilm','Ttn','Lmt','Cpx',...
    'Fe_Amph','Epd_Zo','Apt','Tml','Zrc','Calg','AlOx','Udf','Total'};
mins_abv = {'Qtz','Ksp','Plg','Mscv','Bti','Kao','Chl','Mg Clay','Mg Si',...
    'Ill-Smec','Cal','Dol','Fe-Ox/Sid','Prt','Gyp/Anh','Hlt','Rtl/Ilm',...
    'Ttn','Lmt','Cpx','Fe Amph','Epd/Zo','Apt','Tml','Zrc','AlOx','Tot'};
minerals = var_BMMP_F(1,7:end-1);
mineralsA2 = minerals([3,4,5,6,8,9,10,11,13,14,17,18,20,21,26,27]);
mineralsA4 = [minerals([1,2,3]),'An>25',minerals([7,8]),'Biotite high',...
    minerals(12),'Chlorite',minerals([15,16]),'Illite/smectite',...
    minerals(19),'Dolomite',minerals([22:25]),'Rutile/ilmenite',...
    minerals([28:35,37])];
mineralsA = {'Quartz','K Feldspar','Plagioclase','Muscovite','Biotite','Kaolinite',...
    'Chlorite','Mg Clay','Mg Silicate','Illite-Smectite','Calcite','Dolomite',...
    'Fe-Oxide-Siderite','Pyrite','Gypsum-Anhydrite','Halite','Rutile-Ilmenite',...
    'Titanite','Laumontite','Clinopyroxene','Fe Amphibole','Epidote-Zoisite',...
    'Apatite','Tourmaline','Zircon','Aluminum Oxide'};

elements = varMD(2,6:end)
elements2 = {'Si','Al','Ba','Be','Ca','Cr','Cu','Fe','K','Li','Mg','Mn','Ni',...
    'P','Sc','Ti','V','Zn','As','Cd','Ce','Co','Cs','Dy','Er','Eu','Ga','Gd',...
    'Hf','Ho','La','Lu','Mo','Nb','Nd','Pb','Pr','Rb','Sb','Sm','Tb','Th',...
    'Tm','U','Y','Yb','Zr'}

labels = var_BMMP_F(2:end,1);
ccat = [1,2,[3,4,5,6],7,[8,9,10,11],12,[13,14],15,16,[17,18],19,[20,21],22,23,24,25,[26,27],28,29,30,31,32,33,34,35,37,39];

%% data matrix
% numA = [58,69,80,91]
% numA = [14,25,36,47]
% for j = 1:4

%     dat = dat_SMPG_F(4:4:end,:)
    dat = dat_BMMP_F(keepLab-1,:)

    datMtx = dat   %%%in this grouping you are simply adding up all similar minerals
    abv = 0
    if abv == 1
    mineralOut = mineralsA
    datMtxMin = [dat(:,1:2),sum(dat(:,3:6),2),dat(:,7),sum(dat(:,8:11),2),...
        dat(:,12),sum(dat(:,13:14),2),dat(:,15:16),sum(dat(:,17:18),2),...
        dat(:,19),sum(dat(:,20:21),2),dat(:,22:25),sum(dat(:,26:27),2),...
        dat(:,28:35),dat(:,37)];
    end
    
    abv2 = 0    %%%these are the minerals that you were grouping up, but on their own (biotite, Plag, etc.)
    if abv2 == 1
    mineralOut = mineralsA2
    datMtxMin = dat(:,[3,4,5,6,8,9,10,11,13,14,17,18,20,21,26,27]);
    end
    
    abv3 = 1     %%%heres the full suite without the annoying "total", "udf" "calgon" etc
    if abv3 == 1
        mins = mins([1:35,37])
        mineralOut = minerals([1:35,37])
        datMtxMin = dat(:,[1:35,37]);
    end
    
    abv4 = 0   %%%here you lump like abv 3, but separate out plag from albite and biotite low from int+high
    if abv4 == 1
        mineralOut = mineralsA4;
        datMtxMin = [dat(:,1:3),sum(dat(:,4:6),2),dat(:,7:8),sum(dat(:,9:11),2),...
        dat(:,12),sum(dat(:,13:14),2),dat(:,15:16),sum(dat(:,17:18),2),...
        dat(:,19),sum(dat(:,20:21),2),dat(:,22:25),sum(dat(:,26:27),2),...
        dat(:,28:35),dat(:,37)];
    end
    
    abvChem = 1
    
    if abvChem == 1
        elems = elements2
        j = 1
        el1 = []
        for i = 1:length(elements2)
            el1(i) = find(strcmp(elements,elements2(i)))
        end
        datMtxElm = datMD(:,el1)
    end
    
            
    srtin = 0
    if srtin == 1
        [s,ia] = sort(sum(datMtx,1),'descend')
        datMtx = datMtx(:,ia)
        minerals = minerals(ia)
    end

    return
    %% in this section you can plot all of the minerals together
    % close all
    % f1 = figure
    % for i = 1:length(labels)
    %     
    %     hold on
    %     if p_r(i)
    %         mc = [1 0 0]
    %     elseif ms_r(i)
    %         mc = [0.7 0.7 0.7]
    %     end
    %     p(i) = plot((datMtx(i,:)),'ko','markerfacecolor',mc);
    %     
    % end
    % 
    % XX = 1:length(minerals)
    % set(gca(),'XTick',XX, 'XTickLabel',minerals)
    % rotateXLabels(gca(),60)
    % 
    % xlim([0 length(minerals)+1])
    % ylabel('mass %')
    % grid on
    % legend([p(1) p(2)],{'Plutonic','Metased.'},'location','northeast')
    % savePDFfunction(f1,'QS_fragRock_massP')

    %% in this section you can make 4*2 subplots for each mineral, and do some 
    %%% stats while you're at it
    close all
%%%do you want to include the gabbro?
    no16 = logical(1-strcmp(lab,'16b'))
    yMtx = datMtxElm(no16,:)
    
    groupA_p = []
    rtA_p = []
    gtA_p = []
    group2_p = []
    for i = 1:length(elems)

        close all
    %     a = 1
    %     ii = i+(a-1);
         y = yMtx(:,i);
        [p,anovatab,stats] = anova1(y,group(no16),'off');
        [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1);
        groupA_p(:,i) = c(:,6);
        title(elems(i));

%         %%%Metaseds surge versus non surge
%         [p,anovatab,stats] = anova1(y(strcmp('MS',rt)),gt(strcmp('MS',rt)),'off');
%         [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1);
%         group2_p(3,i) = c(:,6);
%         title(minerals(i));
%         %%%plutonic surge versus non surge
%         [p,anovatab,stats] = anova1(y(strcmp('P',rt)),gt(strcmp('P',rt)),'off');
%         [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1);
%         group2_p(2,i) = c(:,6);
%         title(minerals(i));
%         %%%surge all rock
%         [p,anovatab,stats] = anova1(y(strcmp('S',gt)),rt(strcmp('S',gt)),'off');
%         [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1);
%         group2_p(1,i) = c(:,6);
%         title(minerals(i));
%         %%%non-surge all rock
%         [p,anovatab,stats] = anova1(y(strcmp('NS',gt)),rt(strcmp('NS',gt)),'off');
%         [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1);
%         group2_p(4,i) = c(:,6);
%         title(minerals(i));

        [p,anovatab,stats] = anova1(y,rt(no16),'off');
        [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1);
        rtA_p(:,i) = c(:,6);
        title(elems(i));

%         [p,anovatab,stats] = anova1(y,gt,'off');
%         [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1);
%         gtA_p(:,i) = c(:,6);
%         title(minerals(i));

    end
    %% write out the p-values
%     num1 = numA(j)
    num1 = 3
    xlswrite('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\p_Out.xlsx',groupA_p(2:5,:),2,['B' num2str(num1) ':AB' num2str(num1+3)])
    xlswrite('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\p_Out.xlsx',group2_p,2,['B' num2str(num1+4) ':AB' num2str(num1+7)])
    xlswrite('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\p_Out.xlsx',rtA_p,2,['B' num2str(num1+8) ':AB' num2str(num1+8)])
    xlswrite('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\p_Out.xlsx',gtA_p,2,['B' num2str(num1+9) ':AB' num2str(num1+9)])

% end
%% make some pairwise plots
[var,dat] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\p_Out.xlsx')

% str = 'K Feldspar'
% el = find(strcmp(minerals,str))


vlt = var(1,:)<0.1

% for i = 1:length(vlt)
for i = 1:length(elems)
    
%     if vlt(i)

    %%%two plots of rock type
    close all
    f1 = figure
    y = datMtxElm(:,i)
    str = elems(i)
    hold on
    ms = 10
    fs = 16

    p1 = plot(ones(1,sum(p_r.*st))*4,y(logical(p_r.*st)),'ko','markerfacecolor','r','markersize',ms)
    p2 = plot(ones(1,sum(p_r.*nst))*3,y(logical(p_r.*nst)),'k^','markerfacecolor','r','markersize',ms)
    p3 = plot(ones(1,sum(ms_r.*st))*2,y(logical(ms_r.*st)),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms)
    p4 = plot(ones(1,sum(ms_r.*nst))*1,y(logical(ms_r.*nst)),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms)
    grid on
%     text(p_r+1+0.1,y,lab,'fontsize',fs-2)
    text(groupLab+0.1,y,lab)
    set(gca,'XTick',1:4,'XTickLabel',{'MSNS','MSS','PNS','PS'})
    xlim([0 5])
    
    xlms = get(gca,'xlim')
    ylms = get(gca,'ylim')
    xt = diff(xlms)*0.1
    yt = diff(ylms)*0.9
%     text(xt,yt,['p = ' num2str(var(1,i),1)],'fontsize',fs)
%     legend([p1 p2],{'surge','non-surge'},'location','northeast')
    title(str)
    ylabel('mass %')
    set(gca,'fontsize',fs)
    savePDFfunction(f1,[str{1}])
%     keyboard

%     end
end


%% here you can just make a few colormaps from the correlations to see who's who in the zoo
close all
datIN = [datMtxElm]
datM = datIN(ms_r,:)
cdm = corr(datM)
limC = -0.6
nx = floor((limC+1)*32)
limC2n = -1+(nx/32)
limC2p = -limC2n
cdm(logical(((cdm<limC2p)+(cdm>-limC2n))==2))=0
labels = [elements]
fx = figure
pcolor(cdm)
colorbar
cm = colormap(jet)
lcm = length(cm(nx:end-nx,:))
cm(nx:end-nx,:)=ones([lcm,3])
colormap(cm)
caxis([-1 1])
set(gca,'xtick',1.5:1:length(labels),'xticklabel',labels)
set(gca,'ytick',1.5:1:length(labels),'yticklabel',labels)
rotateXLabels( gca(), 90, labels)
set(gca,'fontsize',12)
% savePDFfunction(fx, 'Min_Cor_P_ms_r')

%% here you can loop through the nested loop to visually inspect correlations
ms = 10
for i = 1:length(labels)-1
    for j = i+1:length(labels)
       
        x = datIN(:,i);
        y = datIN(:,j);
        close all
        f1 = figure;
%         set(gcf,'visible','off')
        hold on
        p1 = plot(x(logical(p_r.*st)),y(logical(p_r.*st)),'ko','markerfacecolor','r','markersize',ms);
        p2 = plot(x(logical(p_r.*nst)),y(logical(p_r.*nst)),'k^','markerfacecolor','r','markersize',ms);
        p3 = plot(x(logical(ms_r.*st)),y(logical(ms_r.*st)),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
        p4 = plot(x(logical(ms_r.*nst)),y(logical(ms_r.*nst)),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
        p4 = plot(x(logical(Gb_r.*nst)),y(logical(Gb_r.*nst)),'k^','markerfacecolor',[0 0 0],'markersize',ms);
        
        xlm = get(gca,'xlim');
        ylm = get(gca,'ylim');
        
        cMS = corr(x(ms_r),y(ms_r));
        pMS = polyfit(x(ms_r),y(ms_r),1);
        yMS = pMS(1).*xlm + pMS(2);
        plot(xlm,yMS,'k')
        
        cP = corr(x(p_r),y(p_r));
        pP = polyfit(x(p_r),y(p_r),1);
        yP = pP(1).*xlm + pP(2);
        plot(xlm,yP,'color','r');
        
        text(0.8, 0.9,['r_{MS} = ' num2str(cMS,2) ],'units','normalized','fontsize',18);
        text(0.8, 0.8,['r_{P} = ' num2str(cP,2) ],'units','normalized','color','r','fontsize',18);
        
        xlim(xlm);
        ylim(ylm);
        
        xlabel(labels(i));
        ylabel(labels(j));
        grid on
        set(gca,'fontsize',18);
%         savePDFfunction(f1,[labels{i} 'VS' labels{j}])
        keyboard

    end
end



%% x-y scatter by string compare
close all
ms = 10

str1 = 'Illite-Smectite'
str2 = 'Muscovite'

i1 = find(strcmp(mineralOut,str1))
i2 = find(strcmp(mineralOut,str2))
%Illite & illite-smectite
% i1 = find(strcmp(mineralOut,'Muscovite'))
% i2 = find(strcmp(mineralOut,'Illite-Smectite'))
% i3 = find(strcmp(mineralOut,'Rutile-Ilmenite'))
% i4 = find(strcmp(mineralOut,'K Feldspar'))
% j1 = find(strcmp(mineralOut,'Laumontite'))
% j2 = find(strcmp(mineralOut,'Plagioclase'))
% j3 = find(strcmp(mineralOut,'Epidote-Zoisite'))
%         x = sum(X(:,[i1,i2,i3,i4]),2);
%         y = sum(X(:,[j1,j2,j3]),2);
%         close all
        x = datMtxMin(:,i1)
        y = datMtxMin(:,i2)
        f1 = figure;
%         set(gcf,'visible','off')
        hold on
%         p1 = plot(x(logical(p_r.*st)),y(logical(p_r.*st)),'ko','markerfacecolor','r','markersize',ms);
%         p2 = plot(x(logical(p_r.*nst)),y(logical(p_r.*nst)),'k^','markerfacecolor','r','markersize',ms);
        p3 = plot(x(logical(ms_r.*st)),y(logical(ms_r.*st)),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
        p4 = plot(x(logical(ms_r.*nst)),y(logical(ms_r.*nst)),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
%         p4 = plot(x(logical(Gb_r.*nst)),y(logical(Gb_r.*nst)),'k^','markerfacecolor',[0 0 0],'markersize',ms);
%       
        text(x(ms_r,:),y(ms_r,:),lab(ms_r))
        text(x(p_r,:),y(p_r,:),lab(p_r))
        text(x(Gb_r,:),y(Gb_r,:),lab(Gb_r))
        xlm = get(gca,'xlim');
        ylm = get(gca,'ylim');
        
        cMS = corr(x(ms_r),y(ms_r));
        pMS = polyfit(x(ms_r),y(ms_r),1);
        yMS = pMS(1).*xlm + pMS(2);
        plot(xlm,yMS,'k')
        
%         cP = corr(x(p_r),y(p_r));
%         pP = polyfit(x(p_r),y(p_r),1);
%         yP = pP(1).*xlm + pP(2);
%         plot(xlm,yP,'color','r');
        
        text(0.8, 0.9,['r_{MS} = ' num2str(cMS,2) ],'units','normalized','fontsize',18);
%         text(0.8, 0.8,['r_{P} = ' num2str(cP,2) ],'units','normalized','color','r','fontsize',18);
        
        xlim(xlm);
        ylim(ylm);
        
        xlabel(mineralOut(i1));
        ylabel(mineralOut(i2));
        grid on
        set(gca,'fontsize',18);
%         savePDFfunction(f1,[labels{i} 'VS' labels{j}])
      


%% PCA
%%%do some PCA on the datas, first stargin with the bulk mineralogy. This seems like a slight waste of time right now. So moving on.

% samps = logical(ones(1,length(lab)))
samps = find(ms_r)
close all
X = [datMtxMin(samps,:),datMtxElm(samps,:)]
z = X-mean(X)
z = z./(repmat(std(X),[size(X,1),1]))
X = z

vars = [mineralOut,elems]

varCov = cov(X);
[ve,va] = eig(varCov);

Sr = X*ve;
Sr = Sr(:,[end:-1:1]);

ve = ve(:,[end:-1:1]);
ve = -ve;
va = sum(va);
va = va([end:-1:1]);

tv = trace(varCov);

pv = cumsum((va/tv)*100);

%%%plot the percent variation of each component
fb = figure
plot(pv,'-o')

%%% plot the loads
fPC1 = figure(10)
subplot(2,1,1)
a = 1
[ves,els] = sort(ve(:,a)')
stem(ves)
text([1:length(vars)],ves,vars(els))
ylabel('Loads on PC 1')
set(gca,'fontsize',14)

fPC2 = figure(11)
subplot(2,1,1)
a = 2
[ves,els] = sort(ve(:,a)')
stem(ves)
text([1:length(vars)],ves,vars(els))
title('PC 2 loads')
set(gca,'fontsize',14)

fPC3 = figure(12)
subplot(2,1,1)
a = 3
stem(ves)
text([1:length(vars)],ves,vars(els))
title('PC 3 loads')
set(gca,'fontsize',14)

%%%plot the scores against eachother
f1 = figure
ms = 10
x = Sr(:,1)
y = Sr(:,2)
hold on
%%%here's how to plot just the MS rocks
% stM = find(ms_r.*st)
% nstM = find(ms_r.*nst)
% stM = logical(st.*ms_r)
% nstM = logical(nst.*ms_r)

stM = logical(st(ms_r))
nstM = logical(nst(ms_r))
p3 = plot(x(stM),y(stM),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
p4 = plot(x(nstM),y(nstM),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(x,y,lab(ms_r))
% p1 = plot(x(logical(p_r.*st)),y(logical(p_r.*st)),'ko','markerfacecolor','r','markersize',ms);
% p2 = plot(x(logical(p_r.*nst)),y(logical(p_r.*nst)),'k^','markerfacecolor','r','markersize',ms);
% p3 = plot(x(logical(ms_r.*st)),y(logical(ms_r.*st)),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
% p4 = plot(x(logical(ms_r.*nst)),y(logical(ms_r.*nst)),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
% p5 = plot(x(logical(Gb_r.*nst)),y(logical(Gb_r.*nst)),'k^','markerfacecolor',[0 0 0],'markersize',ms);
% xlm = get(gca,'xlim')
% xos = (xlm(2)-xlm(1))/40
% text(x+xos,y,lab)
% grid on

ylabel('PC 2')
xlabel('PC 1')
% legend([p1 p2 p3 p4 p5],{'MXS','MXNS','MSS','MSNS','GNS'},'location','southwest')
set(gca,'fontsize',18)
% savePDFfunction(f1,'PCMinsEls_allSamps_no16')


figure(10)
x = Sr(:,1)
subplot(2,1,2)
hold on
p1 = plot(ones(1,length(x(stM))),x(stM),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(stM)))+0.1,x(stM),lab(logical(st.*ms_r)));
p2 = plot(ones(1,length(x(nstM)))*2,x(nstM),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(nstM)))*2+0.1,x(nstM),lab(logical(nst.*ms_r)));
set(gca,'XTick', 1:2,'XTickLabel',{'MSS','MSNS'})
xlim([0 3])
grid on
ylabel('PC 1')
set(gca,'fontsize',14)
% savePDFfunction(fPC1,'MS_mins_PC1_abv1')

figure(11)
x = Sr(:,2)
subplot(2,1,2)
hold on
p1 = plot(ones(1,length(x(stM))),x(stM),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(stM)))+0.1,x(stM),lab(logical(st.*ms_r)));
p2 = plot(ones(1,length(x(nstM)))*2,x(nstM),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(nstM)))*2+0.1,x(nstM),lab(logical(nst.*ms_r)));
set(gca,'XTick', 1:2,'XTickLabel',{'MSS','MSNS'})
xlim([0 3])
grid on
ylabel('PC 1')
set(gca,'fontsize',14)
% savePDFfunction(fPC,'MS_mins_PC1')

figure(12)
x = Sr(:,3)
subplot(2,1,2)
hold on
p1 = plot(ones(1,length(x(stM))),x(stM),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(stM)))+0.1,x(stM),lab(logical(st.*ms_r)));
p2 = plot(ones(1,length(x(nstM)))*2,x(nstM),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(nstM)))*2+0.1,x(nstM),lab(logical(nst.*ms_r)));
set(gca,'XTick', 1:2,'XTickLabel',{'MSS','MSNS'})
xlim([0 3])
grid on
ylabel('PC 1')
set(gca,'fontsize',14)
% savePDFfunction(fPC,'MS_mins_PC1')
%% Factor analysis
clc
close all

textVars = mineralOut
% samps = ones(1,length(lab))
samps = find(ms_r)
m = 2
x = datMtxMin(samps,:)
z = x-mean(x)
z = z./(repmat(std(x),[size(x,1),1]))
Sigma = cov(z)
sigma_ii = diag((diag(Sigma)))
[ve,va] = eig(Sigma)
ve = ve(:,end:-1:1)
va = diag(va)'
va = va(end:-1:1)
L  = ve.*repmat(va,[size(ve,1),1])
L = rotatefactors(L,'method','varimax')
L = L(:,1:m)
spec_v = Sigma - L*L'
f = inv(L'*L)*L'*z'

figure
plot(cumsum(va./sum(va)*100),'-o')

ffl = figure
plot(L(:,1),L(:,2),'o','markeredgecolor',[1 1 1])
text(L(:,1)',L(:,2)',textVars,'fontsize',14)
xlabel('Factor 1')
ylabel('Factor 2')
set(gca,'fontsize',18)
% savePDFfunction(ffl,'factorLoading_Minerals_VMX')


fA1 = figure(10)
subplot(2,1,1)
% suptitle('4 factor model with varimax rotation')
a = 1
[ves,els] = sort(L(:,a)')
stem(ves)
text([1:length(textVars)],ves,textVars(els))
ylabel('Factor 1 loads')
set(gca,'fontsize',14)
% figure
% subplot(3,1,1)
% stem(L(:,1))
% text(1:length(L(:,1)),L(:,1),elems)
% subplot(3,1,2)
% stem(L(:,2))
% text(1:length(L(:,2)),L(:,2),elems)
% subplot(3,1,3)
% stem(L(:,3))
% text(1:length(L(:,3)),L(:,3),elems)

ffs = figure
ms = 10
x = f(1,:);
y = f(2,:);
hold on


%%%here's how to plot just the MS rocks
stM = st(find(ms_r))
nstM = nst(find(ms_r))
p3 = plot(x(stM),y(stM),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
p4 = plot(x(nstM),y(nstM),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);

% p1 = plot(x(logical(p_r.*st)),y(logical(p_r.*st)),'ko','markerfacecolor','r','markersize',ms);
% p2 = plot(x(logical(p_r.*nst)),y(logical(p_r.*nst)),'k^','markerfacecolor','r','markersize',ms);
% p3 = plot(x(logical(ms_r.*st)),y(logical(ms_r.*st)),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
% p4 = plot(x(logical(ms_r.*nst)),y(logical(ms_r.*nst)),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
% p5 = plot(x(logical(Gb_r.*nst)),y(logical(Gb_r.*nst)),'k^','markerfacecolor',[0 0 0],'markersize',ms);
text(x,y,lab(ms_r))
grid on
% legend([p1 p2 p3 p4 p5],{'PS','PNS','MSS','MSNS','GNS'},'location','southwest')
xlabel('Factor 1 scores')
ylabel('Factor 2 scores')
% title('Five factor model with varimax rotation')
set(gca,'fontsize',18)
% savePDFfunction(ffs,'factor1_2_scores_All_5fm')


figure(10)
subplot(2,1,2)
hold on
p1 = plot(ones(1,length(x(stM))),x(stM),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(stM)))+0.1,x(stM),lab(logical(st.*ms_r)));
p2 = plot(ones(1,length(x(nstM)))*2,x(nstM),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
text(ones(1,length(x(nstM)))*2+0.1,x(nstM),lab(logical(nst.*ms_r)));
set(gca,'XTick', 1:2,'XTickLabel',{'MSS','MSNS'})
xlim([0 3])
grid on
ylabel('Factor 1 scores')
set(gca,'fontsize',14)
% savePDFfunction(fA1,'F1S_ALLmins')


return
figure
ms = 10
x = f(1,:);
y = f(3,:);
hold on
%%%here's how to plot just the MS rocks
% stM = st(find(ms_r))
% nstM = nst(find(ms_r))
% p3 = plot(x(stM),y(stM),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
% p4 = plot(x(nstM),y(nstM),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);

p1 = plot(x(logical(p_r.*st)),y(logical(p_r.*st)),'ko','markerfacecolor','r','markersize',ms);
p2 = plot(x(logical(p_r.*nst)),y(logical(p_r.*nst)),'k^','markerfacecolor','r','markersize',ms);
p3 = plot(x(logical(ms_r.*st)),y(logical(ms_r.*st)),'ko','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
p4 = plot(x(logical(ms_r.*nst)),y(logical(ms_r.*nst)),'k^','markerfacecolor',[0.7 0.7 0.7],'markersize',ms);
p5 = plot(x(logical(Gb_r.*nst)),y(logical(Gb_r.*nst)),'k^','markerfacecolor',[0 0 0],'markersize',ms);
text(x,y,lab)
grid on

%% multiple linear regression

ms_r_p10 = logical(1-isnan(P10(ms_r)))

samps1 = find(ms_r)
samps = samps1(ms_r_p10)
P10_msr = P10(ms_r)
P10_g = P10_msr(ms_r_p10)
labels = lab(samps)
vars = mineralOut

listOfmins = {'Quartz', 'K Feldspar', 'Plagioclase', 'Muscovite',...
    'Biotite', 'Calcite','Laumontite'}

for i = 1:length(listOfmins)
    el(i) = find(strcmp(vars,listOfmins(i)))
end


n_s = length(samps)
n_v = length(el)

% i = find(strcmp(vars,'Quartz'))
% ni1 = 1:length(vars);
% ni = ni1(ni1~=i)

% X = zscore(datMtxMin(samps,:));
% % X = (datMtxMin);
% Y = X(:,i);
% Z = X(:,ni);

%%%to link the mineralogy to fracture intensity
Z = (datMtxMin(samps,:))
Y = (P10_g')

betaHat = (inv(Z'*Z))*Z'*Y
yhat = Z*betaHat
err = Y-yhat
s2 = (err'*err)./(n_s-(n_v+1))

covBhat = s2*inv(Z'*Z)
varB = diag(covBhat)'
tt = tinv(1-0.025,n_s-n_v-1)

UL = betaHat+(tt*sqrt(varB))'
LL = betaHat-(tt*sqrt(varB))'
return


R_square = 1- (sum(err.^2)/(sum((Y-mean(Y)).^2)))

[sbh,isbh] = sort((betaHat))

close all
figure
bar(sbh)
set(gca,'xtick',1:length(vars)-1,'xticklabel',vars(isbh))
rotateXLabels(gca,60)


[~,isbh] = sort(abs(betaHat),'descend')
figure
A = isbh(1)
B = isbh(2)
plot3(Z(:,A),Z(:,B),Y,'o')
hold on
plot3(betaHat(A)*Z(:,A),betaHat(B)*Z(:,B),yhat,'o','markerfacecolor','k')
grid on
zlabel('')
text(betaHat(A)*Z(:,A),betaHat(B)*Z(:,B),Y,lab(ms_r))

%% matlab linear regression
%%% YOU NEED FEWER VARIABLES THAN SAMPLES. WHEN YOU RETURN FROM OTTAWA,
%%% USE THE MINERALS THAT ARE MOST SIG IN FA AND PCA TO FIND PROVIDE
%%% SIGNIFICANT BETA VALUES
clc
ms_r_p10 = logical(1-isnan(P10(ms_r)))

close all
% listOfmins = {'Muscovite','Ilmenite','Fe-Illite & illite-smectite',...
%     'Laumontite','Plagioclase An25','Plagioclase An50',...
%     'Biotite (Mg-rich)','Plagioclase An75'}
% vars = {'Msc','Ilm','Fe_Ilt','Laum','Plag_25','Plag_50',...
%     'Bti_Mg','Plag_75','P_10'}
listOfmins = {'Muscovite','Plagioclase','Biotite',...
    'Laumontite','Illite-Smectite','Calcite','K Feldspar'}
vars = {'Msc','plag','Bti','Lmt','Ilt_smct','Cal','Kfs','P_10'}


for i = 1:length(listOfmins)
    el(i) = find(strcmp(mineralOut,listOfmins(i)))
end

samps1 = find(ms_r)
samps = samps1(ms_r_p10)
P10_msr = P10(ms_r)
P10_g = P10_msr(ms_r_p10)
labels = lab(samps)

tbl = array2table([datMtxMin(samps,el),P10_g'],'VariableNames',vars)
lm = fitlm(tbl)
%% logistic regression model
close all
% listOfmins = {'Muscovite','Plagioclase','Biotite',...
%     'Laumontite','Illite-Smectite','Calcite','K Feldspar'}
% vars = {'Msc','plag','Bti','Lmt','Ilt_smct','Cal','Kfs','P_10'}
listOfmins = {'Muscovite','Fe-Illite & illite-smectite',...
    'Illite & illite-smectite','Laumontite','Plagioclase An25','Plagioclase An50',...
    'Biotite (Mg-rich)','Plagioclase An75','Biotite (Low Fe & Mg)','Ilmenite'}
vars = {'Msc','Fe_Ilt','Ilt','Laum','Plag_25','Plag_50',...
    'Bti_Mg','Plag_75','bti_low','ilm','P_10'}
for i = 1:length(listOfmins)
    el(i) = find(strcmp(mineralOut,listOfmins(i)))
end

% Z = datMtxMin(ms_r,el)
Z = zscore(datMtxMin(ms_r,el))

[B,dev,stats] = mnrfit(Z,categorical(gt(ms_r)))

figure()
b = B(2:end)
bar(b)
set(gca,'xtick',1:length(listOfmins),'xticklabel',listOfmins)
LL = stats.beta - 1.96.*stats.se
UL = stats.beta + 1.96.*stats.se
hold on
plot(1:length(listOfmins),LL(2:end),'r+')
plot(1:length(listOfmins),UL(2:end),'r+')
rotateXLabels(gca,60)

clc
stats.p
plt = [stats.p(2:end)<0.01]';
Z_p = Z(:,plt);
Z_1 = [ones(length(Z_p(:,1)),1),Z_p];

e_term = exp(Z_1*[B(1);b(plt)]);
p = e_term./(1+e_term);

[cellstr(num2str(p<1)),gt(ms_r),lab(ms_r)']
listOfmins(plt)















