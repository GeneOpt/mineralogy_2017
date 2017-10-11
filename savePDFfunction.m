function [] = savePDFfunction(fid, filename)

f1 = fid;
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperOrientation','landscape');
set(f1,'PaperUnits','normalized');
set(f1,'PaperPosition', [0 0 1 1]);
print(f1,['D:\Code\Summer_2013_data\mineral_data\figures\2017\sed_mins_4plot\' filename],'-dpdf','-r0')


