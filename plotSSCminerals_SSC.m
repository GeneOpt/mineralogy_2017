clear all
close all

minerals1 = cellstr(['Actinolite  ';'Albite Low  ';'Albite Cal  ';...
'Ankerite-Dol';'Biotite 1M  ';'Calcite     ';'Chlinoclore ';...
'Clinozoisite';'Gypsum      ';'Ilt/Mcsv 1M ';'Ilt/Mscv 2M1';...
'Laumontite  ';'K-Feldspar  ';'Montmorillon';'Quartz      '])


minerals2 =     [{'Actinolite  '};{'Albite Low  '};...
{'Ankerite-Dol'};{'Biotite     '};{'Calcite     '};{'Chlinoclore '};...
{'Gypsum      '};{'Clinozoisite'};{'Illite/Mcsvt'};...
{'Laumontite  '};{'K-Feldspar  '};{'Montmorillon'};{'Quartz      '};...
{'Sphene      '}];


minerals3 =     [{'Quartz      '};{'K-Feldspar  '};...
{'Plagioclase '};{'Actinolite  '};{'Sphene      '};{'Biotite     '};...
{'Chlinoclore '};{'Clinozoisite'};{'Illite/Mcsvt'};...
{'Laumontite  '};{'Montmorill. '};{'Ank/Dolomite'};{'Calicite    '};...
{'Gypsum      '}];

%  
% GL1CL01 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'B5:B19')
% Gl1CL05 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'D5:D19')
% GL1_11Jul13_06 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'F5:F19')
% GL1_15Jul13_12 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'H5:H19')
% Gl1_18Jul13_23 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'J5:J19')
% Gl1_18Jul13_31 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'L5:L19')
% Gl1_23Jul13_42 = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'N5:N19')
% Gl1_23Jul13_42_susp  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'B24:B38')
% Gl1_24Jul13_45  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'D24:D38')
% Gl1_spring  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'Q5:Q19')
% Gl1_13BH10  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'T5:T19')
% Gl1_13BH50  = xlsread('D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx',1,'V5:V19')
% 
% save('modalMinGL1_SSC.mat', 'GL1CL01','GL1_11Jul13_06',...
%     'Gl1CL05','GL1_15Jul13_12','Gl1_18Jul13_23','Gl1_18Jul13_31',...
%     'Gl1_23Jul13_42','Gl1_23Jul13_42_susp','Gl1_24Jul13_45',...
%     'Gl1_spring','Gl1_13BH10','Gl1_13BH50')

load('modalMinGL1_SSC.mat')
load('TSminComp.mat')

actualSamples  =   {'On ice CL01    ','On ice CL05  ',...
'11Jul13 06     ','15Jul13 12     ','18Jul13 23     ',...
'18Jul13 31     ','23Jul13 42     ','23Jul13 42 susp',...
'24Jul13 45     ','Spring         ','13BH10         ',...
'13BH50         ','Thin sect 12/14','Thin sect 12/17'...
'Thin sect 12/25','Ave thin sect  '}

sampleNames  =    {'Supraglacial 1  ','Supraglacial 2  ',...
'SS  July 11 #06 ','SS  July 15 #12 ','SS  July 18 #23 ',...
'SS  July 18 #31 ','SS  July 23 #42 ','SS  July 23 #42b',...
'SS  July 24 #45 ','SS  May 30 #01  ','Borehole 13H10  ',...
'Borehole 13H50  ','Optical    12/14','Optical    12/17'...
'Optical    12/25','Optical avergae '}

GL1SSCminMtx = [GL1CL01,GL1_11Jul13_06,Gl1CL05,...
    GL1_15Jul13_12,Gl1_18Jul13_23,Gl1_18Jul13_31,...
    Gl1_23Jul13_42,Gl1_23Jul13_42_susp,Gl1_24Jul13_45,...
    Gl1_spring,Gl1_13BH10,Gl1_13BH50]



% f1 = figure('units','normalized','outerposition',[0 0 1 1])
f1 = figure
axes('position', [ 0.1 0.2 0.8 0.7])

b = zeros(1,16);
colF = ['r';'r';'r';'k';'g';'b';'b';'y';'c';'m';'g';'g';'k';'k';'k';'k']
col =['r';'r';'b';'b';'b';'b';'w';'b';'b';'m';'g';'g';'k';'k';'k';'k']
mark = ['s';'o';'s';'o';'v';'d';'>';'>';'<';'s';'s';'o';'s';'o';'^';'d']

% for i = 1:2
%     mins = GL1SSCminMtx(:,i)
%     mins(2) = mins(2) + mins(3)
%     mins = [mins(1:2);mins(4:end);0]
%     mins(9) = mins(9)+mins(10)
%     mins = [mins(1:9);mins(11:end)]
%     for g  = 1:length(mins)
%         if mins(g) == 0 
%             mins(g) = NaN
%         end
%     end
%     
%     reMins(1) = mins(13);
%     reMins(2) = mins(11);
%     reMins(3) = mins(2);
%     reMins(4) = mins(1);
%     reMins(5) = mins(14);
%     reMins(6) = mins(4);
%     reMins(7) = mins(6);
%     reMins(8) = mins(8);
%     reMins(9) = mins(9);
%     reMins(10) = mins(10);
%     reMins(11) = mins(12);
%     reMins(12) = mins(3);
%     reMins(13) = mins(5);
%     reMins(14) = mins(7);
%     
%    
%         XX = 1:length(minerals2)
%         b(i) = plot([1:length(minerals3)]+(4/6), reMins, col(i), 'marker',mark(i), 'markersize', 10,...
%             'linestyle', 'none', 'markerfacecolor', colF(i))
%         hold on
%         legend(b(1:i), sampleNames(1:i))
% end

for i = 3:9
    mins = GL1SSCminMtx(:,i)
    mins(2) = mins(2) + mins(3)
    mins = [mins(1:2);mins(4:end);0]
    mins(9) = mins(9)+mins(10)
    mins = [mins(1:9);mins(11:end)]
      for g  = 1:length(mins)
        if mins(g) == 0 
            mins(g) = NaN
        end
      end
      
    reMins(1) = mins(13);
    reMins(2) = mins(11);
    reMins(3) = mins(2);
    reMins(4) = mins(1);
    reMins(5) = mins(14);
    reMins(6) = mins(4);
    reMins(7) = mins(6);
    reMins(8) = mins(8);
    reMins(9) = mins(9);
    reMins(10) = mins(10);
    reMins(11) = mins(12);
    reMins(12) = mins(3);
    reMins(13) = mins(5);
    reMins(14) = mins(7);
    
        XX = 1:length(minerals2)
        b(i) = plot([1:length(minerals3)], reMins, col(i), 'marker',mark(i), 'markersize', 10,...
            'linestyle', 'none', 'markerfacecolor', colF(i))
        hold on
        legend(b(3:i), sampleNames(3:i))
end
% 
% for i = 10
%     mins = GL1SSCminMtx(:,i)
%     mins(2) = mins(2) + mins(3)
%     mins = [mins(1:2);mins(4:end);0]
%     mins(9) = mins(9)+mins(10)
%     mins = [mins(1:9);mins(11:end)]
%      for g  = 1:length(mins)
%         if mins(g) == 0 
%             mins(g) = NaN
%         end
%      end
%      
%     reMins(1) = mins(13);
%     reMins(2) = mins(11);
%     reMins(3) = mins(2);
%     reMins(4) = mins(1);
%     reMins(5) = mins(14);
%     reMins(6) = mins(4);
%     reMins(7) = mins(6);
%     reMins(8) = mins(8);
%     reMins(9) = mins(9);
%     reMins(10) = mins(10);
%     reMins(11) = mins(12);
%     reMins(12) = mins(3);
%     reMins(13) = mins(5);
%     reMins(14) = mins(7);
%     
%         XX = 1:length(minerals2)
%         b(i) = plot([1:length(minerals3)]+(2/6), reMins, col(i), 'marker',mark(i), 'markersize', 10,...
%             'linestyle', 'none', 'markerfacecolor', colF(i))
%         hold on
%         legend(b(1:i), sampleNames(1:i))
% end
% 
% for i = 11:12
%     mins = GL1SSCminMtx(:,i)
%     mins(2) = mins(2) + mins(3)
%     mins = [mins(1:2);mins(4:end);0]
%     mins(9) = mins(9)+mins(10)
%     mins = [mins(1:9);mins(11:end)]
%     
%     for g  = 1:length(mins)
%         if mins(g) == 0 
%             mins(g) = NaN
%         end
%     end
%     
%     reMins(1) = mins(13);
%     reMins(2) = mins(11);
%     reMins(3) = mins(2);
%     reMins(4) = mins(1);
%     reMins(5) = mins(14);
%     reMins(6) = mins(4);
%     reMins(7) = mins(6);
%     reMins(8) = mins(8);
%     reMins(9) = mins(9);
%     reMins(10) = mins(10);
%     reMins(11) = mins(12);
%     reMins(12) = mins(3);
%     reMins(13) = mins(5);
%     reMins(14) = mins(7);
%     
%      
%         XX = 1:length(minerals2)
%         b(i) = plot([1:length(minerals3)]+(3/6), reMins, 'color', [0 0.5 0], 'marker',mark(i), 'markersize', 10,...
%             'linestyle', 'none', 'markerfacecolor', [0 0.5 0])
%         hold on
%         legend(b(1:i), sampleNames(1:i))
% end

% 
% for i = 1:4
%     mins = TSModalComp(i,:)
%     for g  = 1:length(mins)
%         if mins(g) == 0 
%             mins(g) = NaN
%         end
%     end
%     
%     reMins(1) = mins(13);
%     reMins(2) = mins(11);
%     reMins(3) = mins(2);
%     reMins(4) = mins(1);
%     reMins(5) = mins(14);
%     reMins(6) = mins(4);
%     reMins(7) = mins(6);
%     reMins(8) = mins(8);
%     reMins(9) = mins(9);
%     reMins(10) = mins(10);
%     reMins(11) = mins(12);
%     reMins(12) = mins(3);
%     reMins(13) = mins(5);
%     reMins(14) = mins(7);
%     
%         j = i+12
%         b(j) = plot([1:length(minerals3)]+(5/6), reMins, col(j), 'marker',mark(j), 'markersize', 10,...
%             'linestyle', 'none', 'markerfacecolor', colF(j))
%         hold on
%         legend(b(1:j), sampleNames(1:j))
% end

set(gca, 'XTick', XX, 'XTickLabel', '', 'fontsize', 16)
h = get(gca,'xlabel')
set(h, 'Units', 'data')
pos = get(h, 'position')
yy = pos(2)
for i = 1:size(minerals2,1)
    t(i) = text(XX(i), yy-1, minerals3(i,:));
end
set(t,'rotation', 60, 'horizontalalignment', 'right')
set(t,'verticalalignment',...
    'middle', 'fontsize',24)
xlim([0 length(minerals2)+1])
ylim([-2 45])
set(f1,'units','normalized')
set(f1,'outerposition', [0 0 1 1])
yl=ylabel('Modal Percent')
set(yl,'fontsize',24)

grid on

set(f1,'Units','Inches');

pos = get(f1,'Position');

set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(f1,'D:\Documents\Thesis\Template_Latex\figures\minCompSSC','-depsc2','-r0')


