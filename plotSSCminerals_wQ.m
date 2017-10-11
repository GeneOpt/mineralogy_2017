run('plotSSCminerals2.m')
%%%%These names and array and such come from 'plotSSCminerals2.m'
close all

SSCSummer = [
0.384
0.84
0.65
1.028
2.83
2.6
2.04
2.485
2.54
3.075
2.58
1.48
3.11
1.31
1.66
1.29
1.21
1.6
3.0225
2.91
1.6075
1.75
1.18
0.97
1.14]

QSummer = [
0.928869827
1.158869523
1.165146394
1.208649554
2.119589631
2.695493288
2.775761783
2.976174596
1.876175441
2.866420416
3.497704487
4.1395976
4.170917729
2.108439649
2.153083238
1.543379436
2.160773855
2.544249663
2.963653331
2.133169777
1.591502071
1.345886704
1.195357501
1.004089389
0.838504387
]

sampleNames  =    {'Supraglacial 1  ','Supraglacial 2  ',...
'SS  July 11 #06 ','SS  July 15 #12 ','SS  July 18 #23 ',...
'SS  July 18 #31 ','SS  July 23 #42 ','SS  July 24 #45 ',...
'SS  July 23 #42b','SS  May 30 #01  ','Borehole 13H10  ',...
'Borehole 13H50  ','XRD Reitv. 12/14','XRD Reitv. 12/25',...
'Progla. till 1  ','Progla. till 2  '};

mineralsRe =     [{'Quartz      '};{'K-Feldspar  '};...
{'Plagioclase '};{'Actinolite  '};{'Biotite     '};...
{'Chlinoclore '};{'Clinozoisite'};{'Illite/Mcsvt'};...
{'Laumontite  '};{'Montmorill. '};{'Ank/Dolomite'};{'Calicite    '};...
{'Gypsum      '}];

timeClay = [196.375 199.3229 199.7291 204.8125 205.0729 ]
Q = [1.4186, 1.9095, 4.1709, 2.1332, 1.3542]
SSC = [1.028 2.54 3.1 2.9 1.75]
minMtx = reMinsMtx([4 5 6 7 8],:)'/100

for i = 1:length(minMtx(1,:))
    minConcMtx(:,i) = minMtx(:,i)*SSC(i)
end

for i = 1:length(minMtx(1,:))
    minFluxMtx(:,i) = minConcMtx(:,i)*Q(i)
end

SSF = SSC.*Q
% 
% x = Q
% y = SSC.*Q
% 
% plot(x,y,'s')
% 
% hold on
% 
% x = QSummer
% y = SSCSummer.*QSummer
% plot(x,y,'ro')

% hold on
% 
x =Q
y = minMtx(12,:)
plot(x,y,'k*')

% hold on
% 
% y = minConcMtx(5,:)
% plot(x,y,'m*')
% 
% y = minConcMtx(6,:)
% plot(x,y,'g*')



convtoconc = fit(y', x','poly1');

ps=coeffvalues(convtoconc);

AA = ps(1);
B = ps(2);

AA = num2str(AA);
B = num2str(B);

linec = linspace(y(1),y(end),length(y));
linconc = y*ps(1) + ps(2);

rs = rsquare(x, linconc)

rs = num2str(rs)
textcon = (['y = ' AA 'x ' B]);
textrs = (['R-square =' rs]);

hold on 
plot(linconc, y)


% minMtx/100
% 
% x = Q
% y = calcite
% plot(x,y, 'o')
% text(x,y,sampleNames(:,[4 5 6 7 8]))