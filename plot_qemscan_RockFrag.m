run read_in_qemscan

%% plotting specs
mt_pr = 'ko'
mt_ms = 'ko'
mfc_s = [1 0 0]
mfc_ns = [0.7 0.7 0.7]
ms = 2

%% here are the specs for the rockFrags
st = strcmp(var_ss(2:end,3),'S')
nst = strcmp(var_ss(2:end,3),'NS')
gt = var_ss(2:end,3)
p_r = strcmp(var_ss(2:end,4),'P')
ms_r = strcmp(var_ss(2:end,4),'MS')
rt = var_ss(2:end,4)
mx_c = strcmp(var_ss(2:end,5),'MX')
ms_c = strcmp(var_ss(2:end,5),'MS')
group = var_ss(2:end,6)
lab = var_ss(2:end,7)
lab = {'1','2','4','5','6a','6b','7','8','9a','9b','10a','10b',...
    '11a','11b','14a','14b','16a','16b','17','18a','18b','19a','19b',...
    '20a','20b','21'}


%%
mins = {'Qtz','Ksp','An0','An25','An50','An75','Mscv','Bti_low','Bti_Mg','Bti_int',...
    'Bti_Fe','Kao','Chl Fe','Chl Mg','Mg Clay','Mg Si','Ill-Smec','Ill-Smec Fe','Cal',...
    'Dol','Dol Fe','Fe-Ox/Sid','Prt','Gyp/Anh','Hlt','Rtl/Ilm','Ilm','Ttn','Lmt','Cpx',...
    'Fe Amph','Epd/Zo','Apt','Tml','Zrc','Calg','AlOx','Udf','Total'};
mins_abv = {'Qtz','Ksp','Plg','Mscv','Bti','Kao','Chl','Mg Clay','Mg Si',...
    'Ill-Smec','Cal','Dol','Fe-Ox/Sid','Prt','Gyp/Anh','Hlt','Rtl/Ilm',...
    'Ttn','Lmt','Cpx','Fe Amph','Epd/Zo','Apt','Tml','Zrc','AlOx'};
minerals = var_BMMP_F(1,7:end-1);
mineralsA = {'Quartz','K Feldspar','Plagioclase*','Muscovite','Biotite*','Kaolinite',...
    'Chlorite*','Mg Clay','Mg Silicate','Illite-Smectite*','Calcite','Dolomite*',...
    'Fe-Oxide/Siderite','Pyrite','Gypsum/Anhydrite','Halite','Rutile/Ilmentite*',...
    'Titanite','Laumontite','Clinopyrozene','Fe Amphibole','Epidote/Zoisite',...
    'Apatite','Tourmaline','Zircon','Aluminum Oxide'};
labels = var_BMMP_F(2:end,1);
ccat = [1,2,[3,4,5,6],7,[8,9,10,11],12,[13,14],15,16,[17,18],19,[20,21],22,23,24,25,[26,27],28,29,30,31,32,33,34,35,37];

%%
dat = dat_BMMP_F
datMtx = dat
abv = 1
if abv == 1
minerals = mineralsA
datMtx = [dat(:,1:2),sum(dat(:,3:6),2),dat(:,7),sum(dat(:,8:11),2),...
    dat(:,12),sum(dat(:,13:14),2),dat(:,15:16),sum(dat(:,17:18),2),...
    dat(:,19),sum(dat(:,20:21),2),dat(:,22:25),sum(dat(:,26:27),2),...
    dat(:,28:35),dat(:,37)];
end

srtin = 0
if srtin == 1
    [s,ia] = sort(sum(datMtx,1),'descend')
    datMtx = datMtx(:,ia)
    minerals = minerals(ia)
end

%% in this section you can plot all of the minerals together
close all
f1 = figure
for i = 1:length(labels)
    
    hold on
    if p_r(i)
        mc = [1 0 0]
    elseif ms_r(i)
        mc = [0.7 0.7 0.7]
    end
    p(i) = plot((datMtx(i,:)),'ko','markerfacecolor',mc);
    
end

XX = 1:length(minerals)
set(gca(),'XTick',XX, 'XTickLabel',minerals)
rotateXLabels(gca(),60)

xlim([0 length(minerals)+1])
ylabel('mass %')
grid on
legend([p(1) p(2)],{'Plutonic','Metased.'},'location','northeast')
savePDFfunction(f1,'QS_fragRock_massP')

%% in this section you can make 4*2 subplots for each mineral, and do some 
%%% stats while you're at it
close all
for i = 1:length(minerals)
    
%     a = 1
%     ii = i+(a-1)
    y = datMtx(:,i)
%     [p,anovatab,stats] = anova1(y,group)
%     [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1)
%     groupA_p(:,i) = c(:,6)
%     title(minerals(i))
%     keyboard
    
    [p,anovatab,stats] = anova1(y,rt)
    [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1)
    rtA_p(:,i) = c(:,6)
    title(minerals(i))
    keyboard
%     
%     [p,anovatab,stats] = anova1(y,gt)
%     [c,m,h,nms] = multcompare(stats,'ctype','hsd','alpha',0.1)
%     gtA_p(:,i) = c(:,6)
%     title(minerals(i))
%     keyboard

end

%% make some pairwise plots
[var,dat] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\p_Out.xlsx')

% str = 'K Feldspar'
% el = find(strcmp(minerals,str))


vlt = var(1,:)<0.1

for i = 1:length(vlt)
    
    if vlt(i)

    %%%two plots of rock type
    close all
    f1 = figure
    y = datMtx(:,i)
    str = minerals(i)
    hold on
    ms = 10
    fs = 16

    p1 = plot(ones(1,sum(p_r.*st))*2,y(logical(p_r.*st)),mt_pr,'markerfacecolor',mfc_s,'markersize',ms)
    p2 = plot(ones(1,sum(p_r.*nst))*2,y(logical(p_r.*nst)),mt_pr,'markerfacecolor',mfc_ns,'markersize',ms)
    p3 = plot(ones(1,sum(ms_r.*st))*1,y(logical(ms_r.*st)),mt_ms,'markerfacecolor',mfc_s,'markersize',ms)
    p4 = plot(ones(1,sum(ms_r.*nst))*1,y(logical(ms_r.*nst)),mt_ms,'markerfacecolor',mfc_ns,'markersize',ms)
    grid on
    text(p_r+1+0.1,y,lab,'fontsize',fs-2)
    set(gca,'XTick',1:2,'XTickLabel',{'MS','P'})
    xlim([0 3])
    
    xlms = get(gca,'xlim')
    ylms = get(gca,'ylim')
    xt = diff(xlms)*0.1
    yt = diff(ylms)*0.9
    text(xt,yt,['p = ' num2str(var(1,i),1)],'fontsize',fs)
    legend([p1 p2],{'surge','non-surge'},'location','northeast')
    title(str)
    ylabel('mass %')
    set(gca,'fontsize',fs)
    savePDFfunction(f1,['massP_sdbL' str{1}])

    end
end




