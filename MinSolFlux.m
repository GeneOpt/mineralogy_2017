% clear all

load('minMtx2.mat')
load('SSCandQ.mat')

mineralsRe =     [{'Quartz      '};{'K-Feldspar  '};...
{'Plagioclase '};{'Actinolite  '};{'Biotite     '};...
{'Chlinoclore '};{'Clinozoisite'};{'Illite/Mcsvt'};...
{'Laumontite  '};{'Montmorill. '};{'Ank/Dolomite'};{'Calicite    '};...
{'Gypsum      '}];

sampleNames  =    {'Supraglacial 1  ','Supraglacial 2  ',...
'SS  July 11 #06 ','SS  July 15 #12 ','SS  July 18 #23 ',...
'SS  July 18 #31 ','SS  July 23 #42 ','SS  July 24 #45 ',...
'SS  July 23 #42b','SS  May 30 #01  ','Borehole 13H10  ',...
'Borehole 13H50  ','XRD Reitv. 12/14','XRD Reitv. 12/25',...
'Progla. till 1  ','Progla. till 2  '};

minWeight = [60.08, 278.33, 268.62, 853.16, 433.53, 595.22, 454.36, 389.34, 470.44...
             549.07, 206.39, 100.09, 172.17];


can=0.35;
kfs = 1;
mgb=0.5;  %frac
mga=0.3;
mgc=0.7; %frac
nam=0.3; %frac
mgm=1;   %frac
mgd=1;   %frac

%%%[Ca,Mg,K,Na,Si,Al,SO4]
         
chemMins = [[0,0,0,0,1,0,0]',...                           %% Quart
            [0,0,kfs,1-kfs,3,1,0]',...                     %% K-spar    
            [can,0,0,1-can,3-can,1+can,0]',...             %% Plag
            [2,5*mga,0,0,8,0,0]',...                       %% Actinolite  
            [0,3*mgb,1,0,3,1,0]',...                       %% Biotite   
            [0,3*mgc,0,0,3,2,0]',...                       %% Chlinoclore   
            [2,0,0,0,3,3,0]',...                           %% Clinozoisite       
            [0,0,1,0,3,3,0]',...                           %% Illite/Musc   
            [1,0,0,0,4,2,0]',...                           %% Laumontite   
            [0.3*(1-nam),2*mgm,0,0.3*nam,4,2*(1-mgm),0]'...%% Montmorill
            [2*(1-mgd),mgd,0,0,0,0,0]',...                 %% Ank/Dol   
            [1,0,0,0,0,0,0]',...                           %% Calcite
            [1,0,0,0,0,0,1]'];                             %% Gypsum
        

SSCMins = reMinsMtx([3,4,5,6,7,8,10],:)

nanSSCMins = isnan(SSCMins)
for i = 1:length(SSCMins(:,1))
    for j = 1:length(SSCMins(1,:))
       if nanSSCMins(i,j) == 1
           SSCMins(i,j) = 0;
       end
    end
end


SSCMinsAve = sum(SSCMins)/length(SSCMins(:,1))/100

%%% this loops through all sidcharge and SSC values, and multiplies it by
%%% the average SSC by wight% to calculate SSC and flux for mineral by
%%% moles and weight
for i = 1:length(SSC)
    ssc = SSC(i);
    QW = Q(i);
    SSCbyMin(i,:) = ssc*SSCMinsAve;
    QbyMin(i,:) = QW*SSCbyMin(i,:);
    SSCmolesMin(i,:) = SSCbyMin(i,:)./minWeight;
    QmolesMin(i,:) = QbyMin(i,:)./minWeight;
end

%%% This is to calculate the elemental abundance for each SSC measured from
%%% all minerals. the rows are the ion/element type, the columns are the
%%% mineral, while the sheets are the samples. 

for i = 1:length(SSC)
    molMin = [];
    elMinMtx = [];
    for j = 1:length(SSCmolesMin(1,:));
    molMin = SSCmolesMin(i,j);
    elMinMtx(:,j) = molMin*chemMins(:,j);
    end
    elMin3DMtx(:,:,i) = elMinMtx;
    sumElMinMtx(:,i) = sum(elMinMtx')'
end

sumElMinMtx = sumElMinMtx';


%%% now, for every sample, you can pick out the combinations of sediments
%%% that you want to put together to analyzed the flux. Insert the mineral
%%% type/number in the brackets in the j loop
minSum = [8,9]
for k = 1:length(SSC)
    for j = 1:length(minSum)
        sampleMtx = elMin3DMtx(:,:,k);
        l = minSum(j);
        someMins(:,j) = sampleMtx(:,l);
    end
    if length(minSum)>1
        sumSomeMins(:,k) = sum(someMins')';
    else
        sumSomeMins(:,k) = someMins;
    end
end

sumSomeMins = sumSomeMins'*10^6
        

















