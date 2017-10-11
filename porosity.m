%% This is  a script to determine surface area for each mineral based on porosity

clear all
close all

%%%This will generate a Mtx with the columns as mineralsRe and the rows as
%%%sampleNames under the matrux variable reMinsMtx
run('plotSSCminerals2.m')
close
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

sample = 11
reMins = reMinsMtx(sample,1:7)
% reMins(6) = 1.8
reMins = reMins/100

%%% now I need to get the length/radius values in the following order
%minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi']);

run('GL1_1225.m')

reOrd = [1,3,2,8,7,5,9]
%%% they now appear in this order
%%% minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi']);
aveLengthTot = aveLengthTot(reOrd)
%%% qtz,micro,plag,hbl,bt,chl,ep
aveLengthTot(4) = aveLengthTot(1)/2
aveLengthTot(6) = aveLengthTot(1)/2
% aveLengthTot = aveLengthTot/1000

% rPArr = [0.2:0.01:0.8]
rPArr = [0.5]
acrP = 1
density=[2.65,2.56,2.67,3.04,3.09,2.65,3.34] %%% g/cc
%%% this is the mass for a porosity that is ignorant to particle size
% porArr = [0.2:0.01:0.6]/1000
% 
% Vt = 1
% denomV = sum(reMins./density)

% 
% aveLengthTot(5) = rPArr(i)    

%%% this is the wight of each mineral for the desired mass per litre.
rP = 50
acrP = 0.5
Ar = [1,1,1,acrP,rP,1,1]
grams = 1
minWeight = reMins*grams
%%% define the lengths as whatever you want, this is in m
aveLengthTot = ones(1,7)*1e-5 %%%10 micrometers
aveLengthTot = aveLengthTot*100 %% now it is cm
%%% the volume per particle in cm3 is
Vpp = (aveLengthTot.^3)./Ar
%%% the mass per particle in grams is
Mpp = Vpp.*density
%%% the surface area per particle is
Sap = (2+(4./Ar)).*(aveLengthTot.^2)
%%% the number of particles is the total mass / Mpp
Np = minWeight./Mpp
%%% the total surface area is
TSA = Np.*Sap
%%% the surface area per gram, is then
SAPG = TSA./minWeight
sum(minWeight)

%%% Mass of particles

disp('Mass of particles')
disp({'Qtz', 'Micro', 'Plag', 'Hbl', 'Biot', 'Chl', 'Epi'})
disp(minWeight)
disp('cm^2 per gram' )
disp(SAPG)






% 
% por = 0.3
% 
% rP = 20
% Ar = [1,1,1,acrP,rP,1,1]
% 
% Mr = (Vt*(1-por))/denomV
% Vr = reMins./density*Mr
% 
% SA = (2+ (4./Ar)).*((Vr.*Ar).^(2/3))
% 
% %%% This is using the particle radius
% 
% %%% this is the volume per particle
% Vpp = (aveL.^3)./Ar
% %%% so the number of particles is the total volume of each mineral divided
% %%% by the volume per particle
% np = Vr./Vpp
% %%% and the total surface area is number of particles times the SA per
% %%% particle
% SApp = (2+(4./Ar)).*(aveL.^2)
% SAt = SApp.*np
% tSAt = sum(SAt)
% % 
% % S2 = (Vr./aveL).*(Ar+2)
% % S2b(i) = S2(5)/sum(S2)
% % 
% 
% SAb(i) = SAt(5)/tSAt
% SAk(i) = SAt(2)/tSAt
% SAp(i) = SAt(3)/tSAt
% MrArr(i) = Mr
% 
% plot(rPArr,SAb,'k')
% hold on
% plot(rPArr,SAk,'r')
% hold on
% plot(rPArr,SAp,'b')
% xlabel('Average length of biotite')
% ylabel('fraction of mineral area')
% xlim([0.2 0.8])















