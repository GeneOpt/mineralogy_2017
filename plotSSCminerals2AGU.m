%%% everything here is for glacier 1, does not include optical, and reads
%%% from 2013_clayTypes_ordered.xlsx
% 
clear all
close all

minerals = cellstr(['Actinolite  ';'Albite Low  ';'Albite Cal  ';...
'Ankerite-Dol';'Biotite 1M  ';'Calcite     ';'Chlinoclore ';...
'Gypsum      ';'Clinozoisite';'Ilt/Mcsv 1M ';'Ilt/Mscv 2M1';...
'Laumontite  ';'K-Feldspar  ';'Montmorillon';'Quartz      ';...
'Pyrite      ';'Vermiculite ';'Alunite     ';'Talc        ';...
'Clinoptiloli';]);

% minerals = cellstr(['Actinolite  ';'Albite';...
% 'Ankerite-Dol';'Biotite 1M  ';'Calcite     ';'Chlinoclore ';...
% 'Gypsum      ';'Clinozoisite';'Ilt/Mcsv';...
% 'Laumontite  ';'K-Feldspar  ';'Montmorillon';'Quartz      ';...
% 'Pyrite      ';'Vermiculite ';'Alunite     ';'Talc        ';...
% 'Clinoptiloli';]);

mineralsRe =    [{'Plagioclase '};{'Quartz      '};...
{'K-feldspar  '};{'Clinozoisite'};{'Biotite     '};...
{'Actinolite '};{'Clinochlore  '};{'Illite/Mscvt.'};...
{'Laumontite  '};{'Montmorill. '};{'Calcite     '};{'Dolomite    '}];

% GL1CL01 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'B5:B24')
% Gl1CL05 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'D5:D24')
% GL1_11Jul13_06 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'F5:F24')
% GL1_15Jul13_12 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'H5:H24')
% Gl1_18Jul13_23 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'J5:J24')
% Gl1_18Jul13_31 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'L5:L24')
% Gl1_23Jul13_42 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'N5:N24')
% Gl1_23Jul13_42_susp  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'B31:B50')
% Gl1_24Jul13_45  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'D31:D50')
% Gl1_spring  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'Q5:Q24')
% Gl1_13BH10  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'T5:T24')
% Gl1_13BH50  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'V5:V24')
% Gl1_12_14   = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'Q31:Q50')
% Gl1_12_25   = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'R31:R50')
% Gl1_Ti1301  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'AB5:AB24')
% Gl1_Ti1302  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes_ordered.xlsx',1,'AC5:AC24')
% 
% 
% save('modalMinGL1_SSC2.mat', 'GL1CL01','Gl1CL05','GL1_11Jul13_06',...
%     'GL1_15Jul13_12','Gl1_18Jul13_23','Gl1_18Jul13_31',...
%     'Gl1_23Jul13_42','Gl1_23Jul13_42_susp','Gl1_24Jul13_45',...
%     'Gl1_spring','Gl1_13BH10','Gl1_13BH50','Gl1_12_14','Gl1_12_25',...
%     'Gl1_Ti1301','Gl1_Ti1302')

load('modalMinGL1_SSC2.mat')

% sampleNames  =    {'SS  July 11 #06 ','SS  July 15 #12 ','SS  July 18 #23 ',...
% 'SS  July 18 #31 ','SS  July 23 #42 ','SS  July 24 #45 ',...
% 'SS  July 23 #42b','SS  May 30 #01  ','Borehole 13H10  ',...
% 'Borehole 13H50  ','Bedrock 1       ','Bedrock 2       ',...
% 'PG till 1       ','PG till 2       '};

GL1SSCminMtx = [GL1_11Jul13_06,...
    GL1_15Jul13_12,   Gl1_18Jul13_23,  Gl1_18Jul13_31,...
    Gl1_23Jul13_42,  Gl1_24Jul13_45,...
    Gl1_spring,  Gl1_13BH10,  Gl1_13BH50,  Gl1_12_14,...
    Gl1_12_25];

% f1 = figure('units','normalized','outerposition',[0 0 1 1])
f1 = figure;
axes('position', [ 0.1 0.2 0.8 0.7]);
b = zeros(1,16);
colE = ['k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k'];
colF =['b';'b';'b';'b';'b';'b';'b';'b';'b';'k';'k'];
mark = ['s';'s';'s';'s';'s';'s';'s';'^';'^';'s';'s'];

for i = 1:length(GL1SSCminMtx(1,:))
    
    mins = GL1SSCminMtx(:,i);
    mins(2) = mins(2) + mins(3);
    mins = [mins(1:2);mins(4:end)];
    mins(9) = mins(9)+mins(10);
    mins = [mins(1:9);mins(11:14)];
    mins = [mins(1:6);mins(8:13)];
%     mins = [mins(1:6);mins(
    for g  = 1:length(mins)
        if mins(g) == 0 
            mins(g) = NaN;
        end
    end
    reMins(1) = mins(2);
    reMins(2) = mins(12);
    reMins(3) = mins(10);
    reMins(4) = mins(7);
    reMins(5) = mins(4);
    reMins(6) = mins(1);
    reMins(7) = mins(6);
    reMins(8) = mins(8);
    reMins(9) = mins(9);
    reMins(10) = mins(11);
    reMins(11) = mins(5);
    reMins(12) = mins(3);
    colrf = colF(i)

if i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 7 
    colrf = [0.5 .8 1]
end
   
if  i == 8 || i == 9
    colrf = [0.3 .5 1]

end
if i == 10 || i == 11
    colrf = [1 0 0]

end
% if i == 10
%     colrf = [0 0.5 0]
% end

    
        XX = 1:length(mineralsRe);
        b(i) = plot([1:length(mineralsRe)]+(4/6), reMins, 'marker',mark(i), 'markersize', 20,...
            'linestyle', 'none', 'markerfacecolor', colrf,'markeredgecolor',colE(i),'linewidth',2);
        hold on
%         legend(b(1:i), sampleNames(1:i));
        reMinsMtx(i,:) = reMins;

end

fs = 36
set(gca, 'XTick', XX, 'XTickLabel', '', 'fontsize', 16)
h = get(gca,'xlabel');
set(h, 'Units', 'data')
pos = get(h, 'position');
yy = pos(2);
for i = 1:size(mineralsRe,1)
    t(i) = text(XX(i)+0.7, yy-1, mineralsRe(i,:));
end
set(t,'rotation', 60, 'horizontalalignment', 'right')
set(t,'verticalalignment',...
    'middle', 'fontsize',fs)
xlim([1 length(mineralsRe)+1])
ylim([-2 55])
set(gca,'fontsize',fs)
set(f1,'units','normalized')
set(f1,'outerposition', [0 0 1 1])
yl=ylabel('Normalized weight (%)');
set(yl,'fontsize',fs,'position',[-0.01 25])

grid on

% savePDFfunction(f1,'mineralogy')

% set(f1,'Units','Inches');
% 
% pos = get(f1,'Position');
% 
% set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% print(f1,'D:\Code\Summer_2013_data\mineral_data\figures\minComp','-dpdf','-r0')

% save('D:\Code\Summer_2013_data\chem_data\minMtx2.mat','reMinsMtx')
