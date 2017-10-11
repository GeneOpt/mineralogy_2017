%%% Glacier 1 12_25



minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi']);
allMins = 1:length(minerals);

%100 = 10mm

a = [0,18,50,60,85,100];
b = [0,10,20,30,50,70,75,85,95,100];
c = [0,10,20,35,60,65,80,90,100];
d = [0,15,20,35,53,60,65,70,75,90,100];
e = [0,15,52,42,55,70,80,90,100];
f = [0,26,30,40,45,60,70,85,90,100];

a = a/(100);
b = b/(100);
c = c/(100);
d = d/(100);
e = e/(100);
f = f/(100);
%%% this now gives a length in cm, for GWB

aa  = [1,1,3,2,7];
bb  = [3,2,3,3,2,1,2,7,1];
cc  = [2,3,2,3,1,1,2,7];
dd  = [2,1,2,1,1,7,2,1,2,2];
ee  = [3,7,2,1,2,1,2,7];
ff  = [2,7,9,3,2,1,2,1,3];


lines = 6;

aveLengthMtx = zeros(lines,length(allMins));
numLengthMtx = zeros(lines,length(allMins));


for k = 1:lines
 
    if k == 1
        ind = a;
        indMin = aa;
    end
    if k == 2
        ind = b;
        indMin = bb;
    end
    if k == 3
        ind = c;
        indMin = cc;
    end
    if k == 4
        ind = d;
        indMin = dd;
    end
    if k == 5
        ind = e;
        indMin = ee;
    end
    
    if k == 6
        ind = f;
        indMin = ff;
    end

    totLength = zeros(1,length(allMins));
    numLength = zeros(1,length(allMins));

    for i = 1:(length(ind)-1)
        lgth = ind(i+1)-ind(i);
        j = indMin(i);
        totLength(j) = totLength(j)+lgth;
        numLength(j) = numLength(j)+1;
    end

    aveLength = (totLength./numLength)
    in = isnan(aveLength);
    for p = 1:length(allMins)
        if in(p) == 1
            aveLength(p) = 0;
        end
    end
    
   
    aveLengthMtx(k,:) = aveLength;
    numLengthMtx(k,:) = numLength;
    
end

totLength = sum(aveLengthMtx);
allLength = sum(totLength);
modalComp25 = sum(aveLengthMtx)/allLength*100;

minCount = sum(1 - (aveLengthMtx == 0));
aveLengthTot = (totLength./minCount);
alt25 = aveLengthTot;

    in = isnan(aveLengthTot);
    for p = 1:length(aveLengthTot)
        if in(p) == 1;
            alt25(p) = 0
        end
    end





