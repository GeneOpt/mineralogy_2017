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
'Borehole 13H50  ','Optical    12/14','Optical    12/17'...
'Optical    12/25','Optical avergae ','XRD Reitv. 12/14',...
'XRD Reitv. 12/15'};

minerals3 =     [{'Quartz      '};{'K-Feldspar  '};...
{'Plagioclase '};{'Actinolite  '};{'Sphene      '};{'Biotite     '};...
{'Chlinoclore '};{'Clinozoisite'};{'Illite/Mcsvt'};...
{'Laumontite  '};{'Montmorill. '};{'Ank/Dolomite'};{'Calicite    '};...
{'Gypsum      '}];

% original
can=0.4;
kfs = 0.5;
mgb=1;
mga=5;
mgc=2;
nam=0.3;
mgm=1;
mgd=1;


% original
%           [Ca,Mg,K,Na,Si,Al,SO4,0]
%              
chemMins = [[0,0,0,0,1,0,0,2]',...                            %% Quartz
            [0,0,kfs,1-kfs,3,1,0,8]',...                          %% K-spar    
            [can,0,0,1-can,3-can,1+can,0,8]',...              %% Plag
            [2,mga,0,0,8,0,0,14]',...                         %% Actinolite
            [0,0,0,0,0,0,0,0]',...                      `     %% Sphene   
            [0,mgb,1,0,3,1,0,12]',...                         %% Biotite   
            [0,mgc,0,0,3,2,0,10]',...                         %% Chlinoclore   
            [1,0,0,0,1,3,0,13]',...                           %% Clinozoisite       
            [0,0,1,0,3,3,0,12]',...                           %% Illite   
            [1,0,0,0,4,2,0,16]',...                           %% Laumontite   
            [0.3*(1-nam),2*mgm,0,0.3*nam,4,2*(1-mgm),0,14]'...%% Montmorill
            [2*(1-mgd),mgd,0,0,0,0,0,6]',...                  %% Ank/Dol   
            [1,0,0,0,0,0,0,3]',...                            %% Calcite
            [1,0,0,0,0,0,1,6]'];                               %% Gypsum
            
        

%%% Here are the XRD mineral comps from 12/14 and 12/25

reMins2 = [27.7,3.9,50.8,1.8,0,4.5,1.8,6.4,3.1,NaN,NaN,NaN,NaN,NaN;...
           31.9,6.9,50.5,0.7,NaN,4.7,NaN,3.5,1.8,NaN,NaN,NaN,NaN,NaN]   ;  
       
%%% concatonate the matrices 
reMinsMtx = [reMinsMtx;reMins2];
       
        
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
% percMassEl = (MassEl.*numEl)./(MassO);
% percMassO = 1-percMassEl;
% 
% for i = 1:5
%     tot = sum(WRCAMtx(i,:))
%     WRCAMtx(i,:) = WRCAMtx(i,:)/tot*100
% end
% 
% for i = 1:5
%     percMass = percMassEl.*WRCAMtx(i,:);
%     percMassOxy = percMassO.*WRCAMtx(i,:);
%     percMassMtx(i,:) = percMass;
%     percMassMtxO(i,:) = percMassOxy;
% end
% 
% for i = 1:5
% totMassPerO = sum(percMassMtxO(i,:));
% totMassPerOArr(i) = totMassPerO;
% end
% 
% for i = 1:5
% MolMtx(i,:) = percMassMtx(i,:)./MassEl;
% end
% 
% MolOx = totMassPerOArr/15.9994;
% 
% MolMtx = [MolMtx,MolOx'];
% 
% for i = 1:5
% totMols = sum(MolMtx(i,:));
% percMoleMtx(i,:) = (MolMtx(i,:)/totMols)*100;
% end
% 
% xlswrite('D:\Field_data\2012\CHemical comp of GL1 Rocks.XLS',percMoleMtx,1,'B47:M51')
% save('wholeRockChem_molePerOx.mat', 'percMoleMtx')
load('wholeRockChem_molePerOx.mat')


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
moleElXRD_OPT(:,11) = fracIonsMtx(:,8);

ElOrder = [{'Si', ' Al', '     Fe', '  Mn', '  Mg', '  Ca', '  Na',...
    '   K', '  Ti', '   P','   O'}];
% 
% disp('Clay on the Glacier form XRD')
% disp(ElOrder)
% disp(moleElXRD_OPT(1:2,:))
% 
% disp('Clay on the Glacier form Whole Rock Chemical')
% disp(ElOrder)
% disp(percMoleMtx(5,:))


disp('Whole rock chemical analysis of bedrock from GL1')
disp(ElOrder)
disp(percMoleMtx(2,:))
disp(percMoleMtx(4,:))

disp('XRD Reitveld GL1')
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









