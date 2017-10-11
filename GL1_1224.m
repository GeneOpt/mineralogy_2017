%%% Glacier 1 12_24



minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi']);
allMins = 1:length(minerals);

a = [0,23,30,45,48,55,62,80,83,90,100];
b = [0,25,35,50,53,60,65,75,95,100];
c = [0,30,45,48,53,80,90,92,100];
d = [0,10,20,30,65,70,75,85,100];
e = [];
f = [];

aa  = [8,1,8,1,8,2,1,8,1,8];
bb  = [8,8,1,8,2,1,8,1,2];
cc  = [1,8,1,8,8,8,1,8];
dd  = [8,1,8,2,4,2,8,8];
ee  = [];
ff  = [];


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

    aveLength = (totLength./numLength)*0.05;
    in = isnan(aveLength);
    for p = 1:length(allMins)
        if in(p) == 1
            aveLength(p) = 0;
        end
    end
    
   
    aveLengthMtx(k,:) = aveLength;
    numLengthMtx(k,:) = numLength;
    
end

numLengthMtx
aveLengthMtx
totLength = sum(aveLengthMtx)
allLength = sum(totLength)
modalComp24 = sum(aveLengthMtx)/allLength*100

minCount = sum(1 - (aveLengthMtx == 0))
aveLengthTot = (totLength./minCount)
alt24 = aveLengthTot;

    in = isnan(aveLengthTot);
    for p = 1:length(aveLengthTot)
        if in(p) == 1;
            alt24(p) = 0;
        end
    end

alt24


