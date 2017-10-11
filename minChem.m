clear all
clc
load('minMtx.mat')

actualSamples  =   {'On ice CL01    ','On ice CL05  ',...
'11Jul13 06     ','15Jul13 12     ','18Jul13 23     ',...
'18Jul13 31     ','23Jul13 42     ','23Jul13 42 susp',...
'24Jul13 45     ','Spring         ','13BH10         ',...
'13BH50         ','Thin sect 12/14','Thin sect 12/17'...
'Thin sect 12/25','Ave thin sect  ','XRD Reitv 12/14',...
'XRD Reitv 12/15'};

sampleNames  =    {'Supraglacial 1  ','Supraglacial 2  ',...
'SS  July 11 #06 ','SS  July 15 #12 ','SS  July 18 #23 ',...
'SS  July 18 #31 ','SS  July 23 #42 ','SS  July 23 #42b',...
'SS  July 24 #45 ','SS  May 30 #01  ','Borehole 13H10  ',...
'Borehole 13H50  ','Optical    12/14','Optical    12/17',...
'Optical    12/25','Optical avergae ','XRD Reitv. 12/14',...
'XRD Reitv. 12/15'};

minerals3 =     [{'Quartz      '};{'K-Feldspar  '};...
{'Plagioclase '};{'Actinolite  '};{'Sphene      '};{'Biotite     '};...
{'Chlinoclore '};{'Clinozoisite'};{'Illite/Mcsvt'};...
{'Laumontite  '};{'Montmorill. '};{'Ank/Dolomite'};{'Calicite    '};...
{'Gypsum      '}];


reMinsMtx(end,2) = 23.1733;
reMinsMtx(end,3) = 23.5468;

%%% Here are the XRD mineral comps from 12/14 and 12/25

reMins2 = [27.7,3.9,50.8,1.8,0,4.5,1.8,6.4,3.1,NaN,NaN,NaN,NaN,NaN;...
           31.9,6.9,50.5,0.7,NaN,4.7,NaN,3.5,1.8,NaN,NaN,NaN,NaN,NaN]   ;  
       
%%% concatonate the matrices 
reMinsMtx = [reMinsMtx;reMins2];
       
% % 
% % original
% can=0.3;
% mgb=3;
% mga=5;
% mgc=5;
% nam=0.3;
% mgm=1;
% mgd=2;
% 
can=0.45;
kfs = 1;
mgb=1;
mga=2;
mgc=1;
nam=0.3; %frac
mgm=1;   %frac
mgd=1;   %frac


% original
%           [Ca,Mg,K,Na,Si,Al,SO4]
             
chemMins = [[0,0,0,0,1,0,0]',...                           %% Quart
            [0,0,kfs,1-kfs,3,1,0]',...                     %% K-spar    
            [can,0,0,1-can,3-can,1+can,0]',...             %% Plag
            [2,mga,0,0,8,0,0]',...                         %% Actinolite
            [0,0,0,0,0,0,0]',...                      `    %% Sphene   
            [0,mgb,1,0,3,1,0]',...                         %% Biotite   
            [0,mgc,0,0,3,2,0]',...                         %% Chlinoclore   
            [1,0,0,0,3,3,0]',...                           %% Clinozoisite       
            [0,0,1,0,3,3,0]',...                           %% Illite   
            [1,0,0,0,4,2,0]',...                           %% Laumontite   
            [0.3*(1-nam),2*mgm,0,0.3*nam,4,2*(1-mgm),0]'...%% Montmorill
            [2*(1-mgd),mgd,0,0,0,0,0]',...                 %% Ank/Dol   
            [1,0,0,0,0,0,0]',...                           %% Calcite
            [1,0,0,0,0,0,1]'];                              %% Gypsum
            
      
% can=0.3;
% mgb=1;
% mga=2;
% mgc=2;
% nam=0.3; %frac
% mgm=1;   %frac
% mgd=1;   %frac
% 
% chemMins = [[0,0,0,0,1,0,0]',...                           %% Quart
%             [0,0,kfs,1-kfs,3,1,0]',...                           %% K-spar    
%             [can,0,0,1-can,3-can,1+can,0]',...             %% Plag
%             [0,0,0,0,0,0,0]',...                         %% Actinolite
%             [0,0,0,0,0,0,0]',...                      `    %% Sphene   
%             [0,mgb,1,0,3,1,0]',...                         %% Biotite   
%             [0,0,0,0,0,0,0]',...                         %% Chlinoclore   
%             [0,0,0,0,0,0,0]',...                           %% Clinozoisite       
%             [0,0,0,0,0,0,0]',...                           %% Illite   
%             [0,0,0,0,0,0,0]',...                           %% Laumontite   
%             [0,0,0,0,0,0,0]'...                           %% Montmorill
%             [0,0,0,0,0,0,0]',...                         %% Ank/Dol   
%             [0,0,0,0,0,0,0]',...                           %% Calcite
%             [0,0,0,0,0,0,0]'];                             %% Gypsum
%         
for i = 1:length(sampleNames)
    mins = reMinsMtx(i,:);
    for j = 1:length(mins);
        molMin = mins(j);
        mmin = (isnan(molMin));
        for k = 1:length(molMin)
            if mmin(k) == 1
                molMin(k) = 0;
            end
        end
        molChem = chemMins(:,j);
        totChem = molMin*molChem;
        tcin = (isnan(totChem));
        for k = 1:length(tcin)
            if tcin(k) == 1
                totChem(k) = 0;
            end
        end
        totChemMtx(:,j) = totChem;
    end
    AllIons = sum(sum(totChemMtx));
    fracIons = sum(totChemMtx')/AllIons*100;
    fracIonsMtx(i,:) = fracIons;
end

% WRCAMtx = xlsread('D:\Field_data\2012\CHemical comp of GL1 Rocks.XLS',1,'B8:K12') 
% MassO = xlsread('D:\Field_data\2012\CHemical comp of GL1 Rocks.XLS',1,'B13:K13') 
% MassEl = xlsread('D:\Field_data\2012\CHemical comp of GL1 Rocks.XLS',1,'B14:K14') 
% numEl = xlsread('D:\Field_data\2012\CHemical comp of GL1 Rocks.XLS',1,'B15:K15') 
% 
% for i = 1:5
%     tot = sum(WRCAMtx(i,:))
%     WRCAMtx(i,:) = WRCAMtx(i,:)/tot*100
% end
% 
% percMassEl = (MassEl.*numEl)./(MassO)
% for i = 1:5
%     percMass = percMassEl.*WRCAMtx(i,:)
%     percMassMtx(i,:) = percMass
% end
% percMassMtx = percMassMtx/(sum(MassO))*100
% for i = 1:5
%     Moles = percMassMtx(i,:)./MassEl
%     MoleMtx(i,:) = Moles
% end
% 
% totMassEl = sum(MoleMtx')
% 
% for i = 1:5
%     percMoles = MoleMtx(i,:)./totMassEl(i)
%     percMoleMtx(i,:) = percMoles*100 
% end
% 
% xlswrite('D:\Field_data\2012\CHemical comp of GL1 Rocks.XLS',percMoleMtx,1,'B41:L45')
% save('wholeRockChem_molePer.mat', 'percMoleMtx')
load('wholeRockChem_molePer.mat')


% [Ca,Mg,K,Na,Si,Al,SO4]
moleElXRD_OPT(:,1) = fracIonsMtx(:,5);
moleElXRD_OPT(:,2) = fracIonsMtx(:,6);
moleElXRD_OPT(:,3) = zeros(1,18);
moleElXRD_OPT(:,4) = zeros(1,18);
moleElXRD_OPT(:,5) = fracIonsMtx(:,2);
moleElXRD_OPT(:,6) = fracIonsMtx(:,1);
moleElXRD_OPT(:,7) = fracIonsMtx(:,4);
moleElXRD_OPT(:,8) = fracIonsMtx(:,3);
moleElXRD_OPT(:,9) = zeros(1,18);
moleElXRD_OPT(:,10) = zeros(1,18);

ElOrder = [{'Si', ' Al', '     Fe', '  Mn', '  Mg', '  Ca', '  Na',...
    '   K', '  Ti', '   P'}];
% 
% disp('Clay on the Glacier form XRD')
% disp(ElOrder)
% disp(moleElXRD_OPT(1:2,:))
% 
% disp('Clay on the Glacier form Whole Rock Chemical')
% disp(ElOrder)
% disp(percMoleMtx(5,:))

fn = 'D:\Code\Summer_2013_data\mineral_data\min2el.xlsx';

disp('Whole rock chemical analysis of bedrock from GL1')
disp(ElOrder)
disp(percMoleMtx(2,:))
disp(percMoleMtx(4,:))


disp('XRD Reitveld analysis of GL1 Rocks')

disp(moleElXRD_OPT(17,:))
disp(moleElXRD_OPT(18,:))

disp('Bedrock GL1 thin section')
disp(moleElXRD_OPT(13,:))
disp(moleElXRD_OPT(14,:))
disp(moleElXRD_OPT(15,:))
disp(moleElXRD_OPT(16,:))

disp('Suspended sediment in stream from summer from XRD')
disp(ElOrder)
disp(moleElXRD_OPT(3,:))
disp(moleElXRD_OPT(4,:))
disp(moleElXRD_OPT(5,:))
disp(moleElXRD_OPT(6,:))
disp(moleElXRD_OPT(7,:))
disp(moleElXRD_OPT(9,:))

disp('Suspended sediment in stream from summer decanting from XRD')
disp(moleElXRD_OPT(8,:))

disp('Suspended sediment in stream from spring from XRD')
disp(moleElXRD_OPT(10,:))

disp('Suspended sediment in stream from boreholes from XRD')
disp(moleElXRD_OPT(11,:))
disp(moleElXRD_OPT(12,:))

wxlsx = 0
if wxlsx == 1

xlswrite(fn, ElOrder, 'B2:K2')
xlswrite(fn,percMoleMtx(2,:), 'B3:K3')
xlswrite(fn, percMoleMtx(4,:), 'B4:K4')

xlswrite(fn, ElOrder, 'B6:K6')
xlswrite(fn,moleElXRD_OPT(17,:), 'B7:K7')
xlswrite(fn, moleElXRD_OPT(18,:), 'B8:K8')

xlswrite(fn, ElOrder, 'B10:K10')
xlswrite(fn,moleElXRD_OPT(13,:), 'B11:K11')
xlswrite(fn, moleElXRD_OPT(14,:), 'B12:K12')
xlswrite(fn,moleElXRD_OPT(15,:), 'B13:K13')
xlswrite(fn, moleElXRD_OPT(16,:), 'B14:K14')

xlswrite(fn, ElOrder, 'B16:K16')
xlswrite(fn,moleElXRD_OPT(3,:), 'B17:K17')
xlswrite(fn, moleElXRD_OPT(4,:), 'B18:K18')
xlswrite(fn,moleElXRD_OPT(5,:), 'B19:K19')
xlswrite(fn, moleElXRD_OPT(6,:), 'B20:K20')
xlswrite(fn,moleElXRD_OPT(5,:), 'B21:K21')
xlswrite(fn, moleElXRD_OPT(6,:), 'B22:K22')

xlswrite(fn, ElOrder, 'B24:K24')
xlswrite(fn,moleElXRD_OPT(8,:), 'B25:K25')

xlswrite(fn, ElOrder, 'B27:K27')
xlswrite(fn,moleElXRD_OPT(10,:), 'B28:K28')

xlswrite(fn, ElOrder, 'B30:K30')
xlswrite(fn,moleElXRD_OPT(11,:), 'B31:K31')
xlswrite(fn,moleElXRD_OPT(12,:), 'B32:K32')

end














